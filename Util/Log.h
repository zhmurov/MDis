/*
 * Log.h
 *
 *  Created on: Apr 7, 2011
 *      Author: serxa
 */


#pragma once

#include <ostream>
#include <sstream>
#include <cassert>
#include <ctime>
#include <iomanip>

//#define DEBUG
#ifdef DEBUG
# define DPRINTF(format, ...) printf(format, ##__VA_ARGS__)
#else
# define DPRINTF(format, ...) // Words unspoken seal the lie...
#endif


//! Interface of a writable log
class ILog {
public:
  virtual ~ILog() {}
  virtual void Write(const char* message) const = 0;
};

//! Wrapper for ILog, that allow you to use ILog as it were a std::ostream
//! It does a single ILog::Write() in destructor
//! \code
//!   LogStream(log) << "hello world!"
//! \endcode
class LogStream : public std::stringstream {
public:
	LogStream(const ILog* log)
		: log_(log), width_(0) {}

	LogStream(const ILog& log)
		: log_(&log), width_(0) {}

	~LogStream() {
		if (log_)
			log_->Write(str().c_str());
	}

	template <class T>
	LogStream& operator << (const T& t) {
		if (log_) {
			if (width_ != 0 && fields_ > 0) {
				static_cast<std::ostream&>(*this) << std::setw(width_);
				--fields_;
			}
			static_cast<std::ostream&>(*this) << t;
		}
		return *this;
	}
	
	//! After calling this method on LogStream, next @fields (default: 4294967294) supplied 
	//! values will be padded with spaces to be at least @width characters long.
	//! Specifically, it was designed to print out data tables with optional comments
	LogStream& table (int width, int fields = -1) {
		width_ = width;
		fields_ = static_cast<unsigned int>(fields);
		return *this;
	}

private:
	const ILog* log_;
	int width_;
	unsigned int fields_;
};

//! \returns a string with current time in format "[%Y-%m-%d %H:%M:%S]"
inline std::string makeTimePrefix() {
	char time_prefix[1024];
	time_t now = time(0);
	tm atm;
	localtime_r(&now, &atm);
	strftime(time_prefix, sizeof(time_prefix), "[%Y-%m-%d %H:%M:%S] ", &atm);
	return std::string(time_prefix);
}
