#ifndef _LOG_VAR_H
#define _LOG_VAR_H

#define ZERO true
#define NOTZERO false
#define POSITIVE true
#define NEGATIVE false
#define LOG_UNDERFLOW 20
#define EQUAL_TOLERANCE 1e-5

#include <iostream>

typedef float lvFloat;

struct lv {
    lvFloat f;
    bool    i   = NOTZERO;
    bool    sng = POSITIVE;
};

class LogVar {

private:
    lv _logVal;

    static lv logAdd    ( lv log_x, lv log_y );
    static lv logMul    ( lv log_x, lv log_y );
    static lv logDiv    ( lv log_x, lv log_y );
    static lv parseFloat( lvFloat f          );

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

    lvFloat logValue() const;
};

LogVar        operator + ( lvFloat a, const LogVar& b         );
LogVar        operator * ( lvFloat a, const LogVar& b         );
LogVar        operator / ( lvFloat a, const LogVar& b         );
std::ostream& operator <<( std::ostream& os, const LogVar& lv );



#endif