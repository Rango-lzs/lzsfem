#pragma once
#include "LogFileStream.h"

#include <stdio.h>
#include <string>

class LogStream;
class LogFileStream;

//-----------------------------------------------------------------------------
//! Class that is used for logging purposes

//! This class can output to different
//! files at the same time.
//! At this time it outputs data to the screen (stdout) and to an external text file.
//! Note that this class is implemented as a singleton, in other words, only one
//! instance can be created.

class Logfile
{
public:
    enum MODE
    {
        LOG_NEVER = 0,
        LOG_FILE = 1,
        LOG_SCREEN,
        LOG_FILE_AND_SCREEN
    };

public:
    //! constructor
    Logfile();

    //! destructor
    ~Logfile();

    //! open a new logfile
    bool open(const char* szfile);

    //! append to existing file
    bool append(const char* szfile);

    //! formatted printing
    void printf(const char* sz, ...);

    //! print a nice box
    void printbox(const char* sztitle, const char* sz, ...);

    //! set the loggin mode
    MODE SetMode(MODE mode);

    //! get the loggin mode
    MODE GetMode();

    //! flush the logfile
    void flush();

    //! close the logfile
    void close();

    //! return the file name
    const std::string& FileName()
    {
        return m_fp->GetFileName();
    }

    //! returns if the logfile is ready to be written to
    bool is_valid()
    {
        return (m_fp != 0);
    }

    // set the log stream
    void SetLogStream(LogStream* ps)
    {
        m_ps = ps;
    }

    // set the file log stream
    void SetFileStream(LogFileStream* fp)
    {
        m_fp = fp;
    }

    // return the file handle
    operator FILE*()
    {
        return (m_fp ? m_fp->GetFileHandle() : 0);
    }

private:
    //! copy constructor is private so that you cannot create it directly
    Logfile(const Logfile& log)
    {
    }

protected:
    LogFileStream* m_fp;  //!< the actual log file

    LogStream* m_ps;      //!< This stream is used to output to the screen

    MODE m_mode;          //!< mode of log file
};
