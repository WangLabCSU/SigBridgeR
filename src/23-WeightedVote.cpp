#include <Rcpp.h>
#include <cstring>

// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

namespace
{

    inline bool is_positive(SEXP x)
    {
        if (x == NA_STRING)
        {
            return false;
        }
        return std::strcmp(CHAR(x), "Positive") == 0;
    }

    struct TieRandom
    {
        static inline int choose()
        {
            return R::unif_rand() < 0.5 ? 0 : 1;
        }
    };

    struct TieFirst
    {
        static inline int choose()
        {
            return 0;
        }
    };

    struct TieLast
    {
        static inline int choose()
        {
            return 1;
        }
    };

    template <typename TiePolicy>
    CharacterVector weighted_vote_impl(
        CharacterMatrix vote_data,
                        NumericVector weights
    )
    {
        const int nr = vote_data.nrow();
        const int nc = vote_data.ncol();
        if (weights.size() != nc)
        {
            stop("Length of `weights` must equal ncol(vote_data).");
        }
        SEXP vote_sexp = vote_data;
        CharacterVector out(nr);
        for (int i = 0; i < nr; ++i)
        {
            double score_pos = 0.0;
            double score_other = 0.0;
            for (int j = 0; j < nc; ++j)
            {
                const double w = weights[j];
                // NA / NaN has no contribution
                if (ISNAN(w) || w == 0.0)
                {
                    continue;
                }
                SEXP v = STRING_ELT(vote_sexp, i + nr * j);
                // NA vote
                if (v == NA_STRING)
                {
                    continue;
                }
                if (is_positive(v))
                {
                    score_pos += w;
                }
                else
                {
                    // Negative / Neutral / Other / -> Other
                    score_other += w;
                }
            }
            int winner;
            if (score_pos > score_other)
            {
                winner = 0;
            }
            else if (score_other > score_pos)
            {
                winner = 1;
            }
            else
            {
                winner = TiePolicy::choose();
            }
            out[i] = winner == 0 ? "Positive" : "Other";
        }
        return out;
    }

} // namespace


// [[Rcpp::export()]]
CharacterVector weighted_vote_cpp(
    CharacterMatrix vote_data,
    NumericVector weights,
    int ties_method
)
{
    Rcpp::RNGScope scope;
    switch (ties_method)
    {
        case 0:
            return weighted_vote_impl<TieRandom>(vote_data, weights);
        case 1:
            return weighted_vote_impl<TieFirst>(vote_data, weights);
        case 2:
            return weighted_vote_impl<TieLast>(vote_data, weights);
        default:
            stop("Invalid `ties_method`.");
    }
    return CharacterVector();
}