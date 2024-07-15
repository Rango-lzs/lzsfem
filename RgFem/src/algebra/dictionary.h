/*****************************************************************//**
 * \file   dictionary.h
 * \brief
 *
 * \author Leizs
 * \date   September 2023
 *********************************************************************/
#ifndef dictionr_h
#define dictionr_h

#include "error.h"

#include <string>
#include <iosfwd>

namespace fem 
{
	class DataStream;
	class Pair;

	/**
	 * This class implements a linked list whose entries are Pairs (see below).
	 *
	 * Dictionaries are typically used by degrees of freedom for storing their unknowns.
	 * A dictionary stores its pairs in a linked list form. It knows the first
	 * pair (attribute 'first') of the list. It also knows the last one (attribute
	 * 'last') in order to append an additional pair fast.
	 */
	class FEM_EXPORT Dictionary
	{
	protected:
		/// First pair
		Pair* first;
		/// Last pair
		Pair* last;

	public:
		/// Constructor, creates empty dictionary
		Dictionary() : first(NULL), last(NULL) { }
		/// Destructor
		~Dictionary();

		/// Clears the receiver.
		void clear();
		/**
		 * Adds a new Pair with given keyword and value into receiver.
		 * @param aKey key of new pair
		 * @param value value of new pair
		 * @return New Pair with given keyword and value
		 */
		Pair* add(int aKey, double value);
		/**
		 * Returns the value of the pair which key is aKey.
		 * If requested key doesn't exist, it is created with assigned value 0.
		 * @param aKey Key for pair.
		 * @return Reference to value of pair with given key
		 */
		double& at(int aKey);
		double at(int aKey) const;
		/**
		 * Checks if dictionary includes given key
		 * @param aKey Dictionary key.
		 * @return True if receiver contains pair with given key, otherwise false.
		 */
		bool includes(int aKey) const;
		/// Prints the receiver on screen.
		void printYourself();
		/// Formats itself as string.
		void formatAsString(std::string& str);
		/// Returns number of pairs of receiver.
		int giveSize();

		/**
		 * Saves the receiver contends (state) to given stream.
		 * @exception throws an ContextIOERR exception if error encountered
		 */
		void saveContext(DataStream& stream);
		/**
		 * Restores the receiver contents (state) from given stream.
		 * @exception throws an ContextIOERR exception if error encountered
		 */
		void restoreContext(DataStream& stream);

		friend std::ostream& operator << (std::ostream& out, const Dictionary& r);
	};
} // end namespace fem
#endif // dictionr_h
