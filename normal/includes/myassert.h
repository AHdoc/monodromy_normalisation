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

template<typename T>
void pop_front(std::vector<T>& vec)
{
    myassert(!vec.empty(),"!vec.empty()");
    vec.erase(vec.begin());
}

template<typename T>
void pop_back(std::vector<T>& vec)
{
    myassert(!vec.empty(),"!vec.empty()");
    vec.pop_back();
}

template<typename T>
void pop_front_back(std::vector<T>& vec)
{
    myassert(!vec.empty(),"!vec.empty()");
    vec.erase(vec.begin());
    vec.pop_back();
}

#endif
