#include <Rcpp.h>

extern "C" {
#include <Rinternals.h>
}

using namespace Rcpp;

// [[Rcpp::export]]
SEXP promise_root_expr_cpp(SEXP expr, SEXP env, int max_depth = 100)
{
    if (!Rf_isEnvironment(env))
    {
        stop("env must be an environment");
    }
    if (max_depth < 1)
    {
        max_depth = 100;
    }
    SEXP cur_expr = expr;
    SEXP cur_env  = env;
    for (int i = 0; i < max_depth; ++i)
    {
        // 只有符号才继续追踪，例如 fn -> x -> y -> Foo
        if (TYPEOF(cur_expr) != SYMSXP)
        {
            return cur_expr;
        }
        if (!Rf_isEnvironment(cur_env))
        {
            return cur_expr;
        }
        // 在当前环境中查找这个符号，不继承父环境
        SEXP val = Rf_findVarInFrame(cur_env, cur_expr);
        if (val == R_UnboundValue)
        {
            return cur_expr;
        }
        // 如果不是 promise，说明已经到根了，例如 Foo 绑定到函数对象
        if (TYPEOF(val) != PROMSXP)
        {
            return cur_expr;
        }
        // 取 promise 的原始表达式和对应环境
        SEXP next_expr = PREXPR(val);
        SEXP next_env  = PRENV(val);
        if (next_expr == R_NilValue || next_env == R_NilValue)
        {
            return cur_expr;
        }
        cur_expr = next_expr;
        cur_env  = next_env;
    }
    return cur_expr;
}