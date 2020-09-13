#include <stddef.h>
#include "BAMReader.h"

#include "RcppArmadillo.h"
using namespace Rcpp;

const char BAMReader::bamEOF[BAMReader::bamEOFlength+1] =
		"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00";
const char BAMReader::bamGzipHead[BAMReader::bamGzipHeadLength+1] = 
		"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00";

// Constructor
BAMReader::BAMReader() {
    bufferPos = 0;
    bufferMax = 0;
    IS_EOF = 0;
    IS_FAIL = 0;
    IS_LENGTH = 0;
}


void BAMReader::SetInputHandle(std::istream *in_stream) {
	IN = in_stream;
    // get length of file:
    IN->seekg (0, std::ios_base::end);
    IS_LENGTH = IN->tellg();
    
    // Check EOF bit
    IN->seekg (-bamEOFlength, std::ios_base::end);
    
    char check_eof_buffer[BAMReader::bamEOFlength+1];
    IN->read(check_eof_buffer, bamEOFlength);
         
    if(strncmp(check_eof_buffer, bamEOF, bamEOFlength) == 0) {
        EOF_POS = IS_LENGTH - bamEOFlength;
        Rcout << "EOF detected at position: " << EOF_POS << "\n";
    } else {
        Rcout << "EOF bit not detected\n";
        EOF_POS = 0;
    }
    IN->seekg (0, std::ios_base::beg);    
}

int BAMReader::LoadBuffer() {
    
    // read compressed buffer
    if(IN->tellg() >= EOF_POS) {
        IS_EOF = 1;
        return(0);
    } else if(IN->fail()) {
        IS_FAIL = 1;
        return(0);
    }

    stream_int16 myInt16;
    // discard bam gzip head. TODO: check bamGzipHead and return error if wrong
    IN->ignore(bamGzipHeadLength);
    IN->read(myInt16.c, 2);

    
    IN->read(compressed_buffer, myInt16.i + 1 - 2  - bamGzipHeadLength);

    bufferMax = 65536;
    z_stream zs;
    zs.zalloc = NULL;
    zs.zfree = NULL;
    zs.msg = NULL;
    zs.next_in = (Bytef*)compressed_buffer;
    zs.avail_in = myInt16.i + 1 - 2  - bamGzipHeadLength;
    zs.next_out = (Bytef*)buffer;
    zs.avail_out = bufferMax;

    stream_int32 myInt32;
    memcpy(myInt32.c, &compressed_buffer[myInt16.i + 1 - 8],4);

    // cout << "Expected CRC = " << myInt32.i << '\t';
        
    // decompress buffer
    // ret = uncompress((Bytef *)buffer, &bufferMax, (Bytef *)compressed_buffer, myInt16.i - 2 - 8 - bamGzipHeadLength + 1);
    int ret = inflateInit2(&zs, -15);
    ret = inflate(&zs, Z_FINISH);
    ret = inflateEnd(&zs);
    
    bufferMax -= zs.avail_out;

    // cout << "Error code: " << ret << '\t' << "avail_out: " << zs.avail_out << '\t';
   
    // check CRC
    uint32_t crc = crc32(crc32(0L, NULL, 0L), (Bytef*)buffer, bufferMax);
    // cout << "CRC = " << crc << '\t';

    bufferPos = 0;
    
    return(ret);
}

int BAMReader::read(char * dest, unsigned int len) {
    
    unsigned int remaining_bytes = 0;
    
    // Initialisation if buffer empty
    if(bufferMax == 0) {
        LoadBuffer();        
    }
    
    if (len < bufferMax - bufferPos) {
        memcpy(&dest[0], &buffer[bufferPos], len);
        bufferPos += len;
        return(Z_OK);
    } else {
        remaining_bytes = len - (bufferMax - bufferPos);
        memcpy(&dest[0], &buffer[bufferPos], bufferMax - bufferPos);
        bufferMax = 0;
        bufferPos = 0;
        
        LoadBuffer();
        
        memcpy(&dest[bufferMax - bufferPos], &buffer[bufferPos], remaining_bytes);
        bufferPos += remaining_bytes;
    }
    return(0);
}

int BAMReader::ignore(unsigned int len) {
    
    unsigned int remaining_bytes = 0;
    
    if (len < bufferMax - bufferPos) {
        // memcpy(dest, &buffer[bufferPos], len);
        bufferPos += len;
        return(Z_OK);
    } else {
        remaining_bytes = len - (bufferMax - bufferPos);
        // memcpy(dest, &buffer[bufferPos], bufferMax - bufferPos);
        bufferMax = 0;
        bufferPos = 0;
        
        LoadBuffer();
        
        // memcpy(dest + bufferMax - bufferPos, &buffer[bufferPos], remaining_bytes);
        bufferPos += remaining_bytes;
    }
    return(0);
}

bool BAMReader::eof() {
    if(IS_EOF == 1) {
        return (true);
    } else {
        if(IN->eof()) {
            IS_EOF = 1;
            return (true);
        } else {
            return (false);
        }
    }
}

bool BAMReader::fail() {
    return(IN->fail());
}

streamsize BAMReader::gcount() {
    return(IN->gcount());
}