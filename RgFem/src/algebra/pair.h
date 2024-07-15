/*****************************************************************//**
 * \file   pair.h
 * \brief
 *
 * \author Leizs
 * \date   September 2023
 *********************************************************************/

#ifndef pair_h
#define pair_h

#include "fem_export.h"

#include <cstdio>

namespace oofem {
	/**
	 * This class implements key/value associations - the key and its associated value.
	 * An instance of Pair is used as an entry in a dictionary.
	 * Pair has three components - its key, its value and pointer to the next Pair in the dictionary.
	 *
	 * Tasks:
	 * - Returning its key, or its value, or the next pair ;
	 * - Appending another pair to itself.
	 *
	 */
	class FEM_EXPORT Pair
	{
	private:
		/// Key.
		int key;
		/// Associate value.
		double value;
		/// Pointer to the next Pair.
		Pair* next;

	public:
		/// Constructor - creates the new Pair with given key k and value v.
		Pair(int k, double v) : key(k), value(v), next(NULL) { }
		/// Destructor
		~Pair() { }

		/// Appends a given pair to itself (sets the pointer to next pair to given Pair).
		void append(Pair* p) { next = p; }
		/// Returns the receiver key.
		int giveKey() { return key; }
		/// Returns pointer to the next pair.
		Pair* giveNext() { return next; }
		/// Returns associated value.
		double& giveValue() { return value; }
		/// Prints receiver to screen.
		void printYourself() { printf("   Pair (%d,%f)\n", key, value); }
	};
} // end namespace oofem
#endif // pair_h
