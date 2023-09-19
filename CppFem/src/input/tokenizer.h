/*****************************************************************//**
 * \file   tokenizer.h
 * \brief  
 * 
 * \author Leizs
 * \date   September 2023
 *********************************************************************/

#ifndef tokenizer_h
#define tokenizer_h

#include "fem_export.h"

#include <vector>
#include <string>

namespace fem {
/**
 * String bracket- and quotation-aware string tokenizer.
 * This class splits given record (represented as string) to particular tokens, which are
 * separated by white spaces.
 * Tokenizer recognizes "quoted strings" and structured tokens that are
 * bounded by '{}, $$' pairs, can be nested and represent single token.
 */
class FEM_EXPORT Tokenizer
{
private:
    /// Array of tokens
    std :: vector< std :: string >tokens;

public:
    /// Constructor. Creates tokenizer with given character as separator.
    Tokenizer();
    /// Tokenizes given record (string).
    void tokenizeLine(const std :: string &line);
    /// returns the number of tokens.
    int giveNumberOfTokens();
    /// Returns pointer to i-th token.
    const char *giveToken(int i);
    //std::string giveToken(int i);

protected:
    /**
     * Reads next simple token (stops when whitespace character is reached)
     * @param pos Starting position.
     * @param line Record from which token is parsed.
     */
    std :: string readSimpleToken(std :: size_t &pos, const std :: string &line);
    /**
     * Reads next token (stops when separator is reached)
     * @param pos Starting position.
     * @param line Record from which token is parsed.
     * @param sep Separator.
     */
    std :: string readToken(std :: size_t &pos, const std :: string &line, char sep);
    /**
     * Reads next structured token (bounded by '{' '}' pairs, possibly nested).
     * @param pos Starting position (should point to a '{').
     * @param line Record from which token is parsed.
     */
    std :: string readStructToken(std :: size_t &pos, const std :: string &line);
    /**
     * Reads next string token (quoted).
     * @param pos Position (index) in token buffer.
     * @param line Record from which token is parsed.
     */
    std :: string readStringToken(std :: size_t &pos, const std :: string &line);
    /**
     * Reads next simple expression token (section identified by starting with '$' and finishing with '$').
     * @param pos Position (index) in token buffer.
     * @param line Record from which token is parsed.
     */
    std :: string readSimpleExpressionToken(std :: size_t &pos, const std :: string &line);
};
} // end namespace fem
#endif // tokenizer_h
