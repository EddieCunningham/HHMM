#ifndef _LOG_VAR_H
#define _LOG_VAR_H

#define ZERO true
#define NOTZERO false
#define UNDERFLOW 500
#define EQUAL_TOLERANCE 1e-5

#include <iostream>

typedef float lvFloat;

struct lv {
    lvFloat f;
    bool    i = NOTZERO;
};

class LogVar {

private:
    lv _logVal;

    static lv _logAdd    ( lv log_x, lv log_y );
    static lv _logMul    ( lv log_x, lv log_y );
    static lv _logDiv    ( lv log_x, lv log_y );
    static lv _parseFloat( lvFloat f          );

public:

    LogVar();
    LogVar( const LogVar& a            );
    LogVar( lv a                       );
    LogVar( lvFloat a                  );
    LogVar( lvFloat a, bool alreadyLog );

    void    operator = ( const LogVar& a );
    void    operator = ( lvFloat a       );

    bool    operator < ( const LogVar& a ) const;
    bool    operator < ( lvFloat a       ) const;

    bool    operator >=( const LogVar& a ) const;
    bool    operator >=( lvFloat a       ) const;

    bool    operator ==( const LogVar& a ) const;
    bool    operator ==( lvFloat a       ) const;

    bool    operator !=( const LogVar& a ) const;
    bool    operator !=( lvFloat a       ) const;

    LogVar  operator + ( const LogVar& a ) const;
    LogVar  operator + ( lvFloat a       ) const;

    LogVar& operator +=( const LogVar& a );
    LogVar& operator +=( lvFloat a       );

    LogVar  operator * ( const LogVar& a ) const;
    LogVar  operator * ( lvFloat a       ) const;

    LogVar& operator *=( const LogVar& a );
    LogVar& operator *=( lvFloat a       );

    LogVar  operator / ( const LogVar& a ) const;
    LogVar  operator / ( lvFloat a       ) const;

    LogVar& operator /=( const LogVar& a );
    LogVar& operator /=( lvFloat a       );

    lvFloat toFloat() const;
};

LogVar        operator + ( lvFloat a, const LogVar& b         );
LogVar        operator * ( lvFloat a, const LogVar& b         );
LogVar        operator / ( lvFloat a, const LogVar& b         );
std::ostream& operator <<( std::ostream& os, const LogVar& lv );



#endif