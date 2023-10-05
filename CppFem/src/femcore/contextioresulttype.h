

#ifndef contextioresulttype_h
#define contextioresulttype_h

namespace fem
{
	enum contextIOResultType {
		CIO_OK = 0,     ///< OK.
		CIO_BADVERSION, ///< Incompatible context file.
		CIO_BADOBJ,     ///< Bad object passed.
		CIO_IOERR       ///< General IO error.
	};
} // end namespace fem
#endif // contextioresulttype_h
