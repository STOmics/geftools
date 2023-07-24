#ifndef GEFTOOLS__FORMAT_H_
#define GEFTOOLS__FORMAT_H_

#include <stdlib.h>

#include <algorithm>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

namespace util {
    class ArgBase {
      public:
        ArgBase() {}
        virtual ~ArgBase() {}
        virtual void Format(std::ostringstream& ss, const std::string& fmt) = 0;
    };

    template <class T>
    class Arg : public ArgBase {
      public:
        Arg(T arg) : m_arg(arg) {}
        virtual ~Arg() {}
        virtual void Format(std::ostringstream& ss, const std::string& fmt) { ss << m_arg; }

      private:
        T m_arg;
    };

    class ArgArray : public std::vector<ArgBase*> {
      public:
        ArgArray() {}
        ~ArgArray() {
            std::for_each(begin(), end(), [](ArgBase* p) { delete p; });
        }
    };

    static void FormatItem(std::ostringstream& ss, const std::string& item, const ArgArray& args) {
        int index = 0;
        int alignment = 0;
        std::string fmt;

        char* endptr = nullptr;
        index = strtol(&item[0], &endptr, 10);
        if (index < 0 || index >= args.size()) {
            return;
        }

        if (*endptr == ',') {
            alignment = strtol(endptr + 1, &endptr, 10);
            if (alignment > 0) {
                ss << std::right << std::setw(alignment);
            } else if (alignment < 0) {
                ss << std::left << std::setw(-alignment);
            }
        }

        if (*endptr == ':') {
            fmt = endptr + 1;
        }

        args[index]->Format(ss, fmt);

        return;
    }

    template <class T>
    static void Transfer(ArgArray& argArray, T t) {
        argArray.push_back(new Arg<T>(t));
    }

    template <class T, typename... Args>
    static void Transfer(ArgArray& argArray, T t, Args&&... args) {
        Transfer(argArray, t);
        Transfer(argArray, args...);
    }

    template <typename... Args>
    std::string Format(const std::string& format, Args&&... args) {
        if (sizeof...(args) == 0) {
            return format;
        }

        ArgArray argArray;
        Transfer(argArray, args...);
        size_t start = 0;
        size_t pos = 0;
        std::ostringstream ss;
        while (true) {
            pos = format.find('{', start);
            if (pos == std::string::npos) {
                ss << format.substr(start);
                break;
            }

            ss << format.substr(start, pos - start);
            if (format[pos + 1] == '{') {
                ss << '{';
                start = pos + 2;
                continue;
            }

            start = pos + 1;
            pos = format.find('}', start);
            if (pos == std::string::npos) {
                ss << format.substr(start - 1);
                break;
            }

            FormatItem(ss, format.substr(start, pos - start), argArray);
            start = pos + 1;
        }

        return ss.str();
    }
}  // namespace util

using LogPrintfFun = std::function<void(const std::string& content)>;
template <typename T>
struct is_char_ptr : std::false_type {};

template <>
struct is_char_ptr<const char*> : std::true_type {};

template <>
struct is_char_ptr<char*> : std::true_type {};

template <typename T>
struct is_null_ptr {
    static bool is(const T& log) { return false; }
};

template <typename T>
struct is_null_ptr<T*> {
    static bool is(const T* log) { return log == nullptr; }
};
void PrintLog(const std::string& content);

class __logwriter {
  public:
    __logwriter(LogPrintfFun printf_fun) : log_printf_fun_(std::move(printf_fun)) {}

    virtual ~__logwriter() {
        if (log_printf_fun_) {
            log_printf_fun_(write_.str());
        }
    }

    template <typename T, typename = typename std::enable_if<!std::is_enum<T>::value>::type>
    inline __logwriter& operator<<(const T& log) {
        if (std::is_pointer<T>::value) {
            if (is_null_ptr<T>::is(log)) {
                write_ << "nullptr";
            } else {
                if (is_char_ptr<T>::value) {
                    write_ << log;
                } else {
                    auto old_flags = write_.flags();
                    write_ << std::hex << log;
                    write_.flags(old_flags);
                }
            }
        } else {
            write_ << log;
        }
        return *this;
    }

    inline __logwriter& operator<<(std::ostream& (*log)(std::ostream&)) {
        write_ << log;
        return *this;
    }

    template <typename T, typename = typename std::enable_if<std::is_enum<T>::value>::type>
    inline __logwriter& operator<<(T v) {
        write_ << static_cast<typename std::underlying_type<T>::type>(v);
        return *this;
    }

  protected:
    std::ostringstream write_;
    LogPrintfFun log_printf_fun_;
};

#define log_info __logwriter(PrintLog)

#endif