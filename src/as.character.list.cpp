#include <Rcpp.h>
#include <sstream>
using namespace Rcpp;

void sexp_write(std::ostringstream& oss, SEXP x);

std::string type_or_class(SEXP x)
{
    SEXP cls = Rf_getAttrib(x, R_ClassSymbol);
    if (!Rf_isNull(cls) && Rf_length(cls) > 0)
    {
        CharacterVector cv(cls);
        std::string out = "object of class ";
        for (int i = 0; i < cv.size(); ++i)
        {
            if (i > 0) out += "/";
            if (cv[i] == NA_STRING)
            {
                out += "NA";
            }
            else
            {
                out += as<std::string>(cv[i]);
            }
        }
        return out;
    }
    switch (TYPEOF(x))
    {
        case NILSXP:  return "NULL";
        case SYMSXP:  return "symbol";
        case LISTSXP: return "pairlist";
        case CLOSXP:  return "function";
        case ENVSXP:  return "environment";
        case PROMSXP: return "promise";
        case LANGSXP: return "call";
        case SPECIALSXP: return "special function";
        case BUILTINSXP: return "builtin function";
        case CHARSXP: return "CHARSXP";
        case LGLSXP: return "logical";
        case INTSXP: return "integer";
        case REALSXP: return "numeric";
        case CPLXSXP: return "complex";
        case STRSXP: return "character";
        case DOTSXP: return "...";
        case ANYSXP: return "any";
        case VECSXP: return "list";
        case EXPRSXP: return "expression";
        case BCODESXP: return "bytecode";
        case EXTPTRSXP: return "external pointer";
        case WEAKREFSXP: return "weak reference";
        case RAWSXP: return "raw";
        case S4SXP: return "S4 object";
        default: return "unknown type";
    }
}

void atomic_write(std::ostringstream& oss, SEXP x)
{
    CharacterVector cv = as<CharacterVector>(x);
    for (int i = 0; i < cv.size(); ++i)
    {
        if (i > 0) oss << ", ";
        if (cv[i] == NA_STRING)
        {
            oss << "NA";
        }
        else
        {
            oss << as<std::string>(cv[i]);
        }
    }
}

void list_write(std::ostringstream& oss, List x)
{
    CharacterVector nms = x.names();
    int n = x.size();
    for (int i = 0; i < n; ++i)
    {
        if (i > 0) oss << ", ";
        if (
            nms.size() == n &&
            nms[i] != NA_STRING &&
            std::string(nms[i]).size() > 0
        )
        {
            oss << as<std::string>(nms[i]);
        }
        else
        {
            oss << i + 1;
        }
        oss << ": ";
        SEXP elem = x[i];
        if (Rf_isNewList(elem))
        {
            oss << "(";
            list_write(oss, as<List>(elem));
            oss << ")";
        }
        else
        {
            sexp_write(oss, elem);
        }
    }
}

void sexp_write(std::ostringstream& oss, SEXP x)
{
    if (Rf_isNull(x))
    {
        oss << "NULL";
        return;
    }
    if (Rf_isNewList(x))
    {
        oss << "(";
        list_write(oss, as<List>(x));
        oss << ")";
        return;
    }
    if (Rf_isVectorAtomic(x))
    {
        atomic_write(oss, x);
        return;
    }
    oss << "<" << type_or_class(x) << ">";
}

// [[Rcpp::export]]
std::string list_to_character_recursive(List x)
{
    std::ostringstream oss;
    list_write(oss, x);
    return oss.str();
}