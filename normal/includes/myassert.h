#ifndef MYASSERT_H
#define MYASSERT_H

#ifndef NDEBUG
#   define myassert(Expr, Msg) \
    __myassert(#Expr, Expr, __FILE__, __LINE__, Msg)
#else
#   define myassert(Expr, Msg) ;
#endif

void __myassert(const char* expr_str, bool expr, const char* file, int line, const char* msg){
	if (!expr){
		std::cerr << "Assert failed:\t" << msg << "\n" << "Expected:\t" << expr_str << "\n" << "Source:\t\t" << file << ", line " << line << "\n";
		abort();
	}
}

#endif
