/**
 * @file backtest_result.hpp
 * @brief Convenience header for BacktestResult.
 *
 * Forwards to backtest_engine.hpp where BacktestResult is defined.
 * This allows the analytics layer to include just the result struct
 * by its logical name without depending on the full engine header.
 */

#ifndef PORTFOLIO_BACKTEST_BACKTEST_RESULT_HPP
#define PORTFOLIO_BACKTEST_BACKTEST_RESULT_HPP

#include "backtest/backtest_engine.hpp"

#endif // PORTFOLIO_BACKTEST_BACKTEST_RESULT_HPP