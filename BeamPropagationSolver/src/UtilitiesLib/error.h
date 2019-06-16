#ifndef EXCEPTION_H
#define EXCEPTION_H

#ifdef DEBUG
	#define Assert(condition,message)    \
	{                                    \
		if(!(condition))                 \
			throw(std::string(message)); \
	}
#else
	#define Assert(condition,message) {}
#endif

#include <mutex>
#include <exception>

class ThreadException {
public:
	ThreadException() : ptr(nullptr) {}

	void rethrow(){
		if(this->ptr)
			std::rethrow_exception(this->ptr);
	}
	void capture_exception() {
		std::unique_lock<std::mutex> guard(this->lock);
		this->ptr = std::current_exception();
	}

private:
	std::exception_ptr ptr;
	std::mutex         lock;
};

#endif
