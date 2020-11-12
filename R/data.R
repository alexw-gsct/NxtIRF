#' Default / recommended filters for NxtIRF
#'
#' This is 
#'  diamonds. The variables are as follows:
#'
#' @format A list of 8 filters:
#'   - Depth filter (requires total splice depth of > 20 (splice + IR depth for IR events))
#'   - Coverage filter IR (90% coverage for IR events)
#'   - Coverage filter AS (binary isoforms tested must represent > 60% of total events across given region)
#'   - Consistency filter MXE / SE (for tandem splice events, each event must be at least 1/4 of total)
#' @md
"default_filter"