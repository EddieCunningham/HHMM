#include "LogVar.h"
#include <cassert>
#include <cmath>

lvFloat log1pExp( lvFloat f ) {
    /* Optimize this later */
    return log1p( exp( f ) );
}

/* -------------------------------------------------------------- */

lv LogVar::_logAdd( lv logX, lv logY ) {

    if( logX.i == ZERO ) { return logY; }
    if( logY.i == ZERO ) { return logX; }

    /* Don't do computations for adding what is practically 0 */
    if( logX.f - logY.f > UNDERFLOW ) { return logX; }
    if( logY.f - logX.f > UNDERFLOW ) { return logY; }

    /* To make sure we'd get underflow instead of overflow */
    /* in the event that |logY.f - logX.f| is large        */
    lv ans;
    if( logX.f > logY.f ) {
        ans.f = logX.f + log1pExp( logY.f - logX.f );
    }
    else {
        ans.f = logY.f + log1pExp( logX.f - logY.f );
    }
    return ans;
}

lv LogVar::_logMul( lv logX, lv logY ) {

    lv ans;
    if( logX.i == ZERO || logY.i == ZERO ) {
        ans.i = ZERO;
    }
    else {
        ans.f = logX.f + logY.f;
    }
    return ans;
}

lv LogVar::_logDiv( lv logX, lv logY ) {

    lv ans;
    if( logX.i == ZERO ) {
        ans.i = ZERO;
    }
    else if( logY.i == ZERO ) {
        /* Can't divide by 0! */
        assert(0);
    }
    else {
        ans.f = logX.f - logY.f;
    }
    return ans;
}

lv LogVar::_parseFloat( lvFloat f ) {

    /* Catch case where f == 0 */
    lv ans;
    if( FP_ZERO == fpclassify( f ) ) {
        ans.i = ZERO;
    }
    else {

        /* Only support positive numbers for the moment */
        assert( f > 0 );

        ans.f = log( f );
    }
    return ans;
}

/* -------------------------------------------------------------- */

LogVar::LogVar() {
    _logVal.i = ZERO;
}

LogVar::LogVar( const LogVar& a ) {
    _logVal = a._logVal;
}

LogVar::LogVar( lv a ) {
    _logVal = a;
}

LogVar::LogVar( lvFloat a ) {
    _logVal = LogVar::_parseFloat( a );
}

LogVar::LogVar( lvFloat a, bool alreadyLog ) {
    if( alreadyLog ) {
        _logVal.i = NOTZERO;
        _logVal.f = log( a );
    }
    else {
        /* Safely parse a float */
        _logVal = LogVar::_parseFloat( a );
    }
}

/* -------------------------------------------------------------- */

void LogVar::operator =( const LogVar& a ) {
    _logVal = a._logVal;
}

void LogVar::operator =( lvFloat a ) {
    _logVal = LogVar( a )._logVal;
}

/* -------------------------------------------------------------- */

bool LogVar::operator <( const LogVar& a ) const {
    if( a._logVal.i == ZERO ) {
        return false;
    }
    if( _logVal.i == ZERO ) {
        return true;
    }
    if( *this == a ) {
        return false;
    }
    return _logVal.f < a._logVal.f;
}

bool LogVar::operator <( lvFloat a ) const {
    LogVar other( a );
    return *this < other;
}

/* -------------------------------------------------------------- */

bool LogVar::operator >=( const LogVar& a ) const {
    return !( *this < a );
}

bool LogVar::operator >=( lvFloat a ) const {
    LogVar other( a );
    return *this >= other;
}

/* -------------------------------------------------------------- */

bool LogVar::operator ==( const LogVar& a ) const {
    if( _logVal.i == ZERO && a._logVal.i == ZERO ) {
        return true;
    }
    if( _logVal.i == NOTZERO && a._logVal.i == NOTZERO ) {
        return std::abs( _logVal.f - a._logVal.f ) < EQUAL_TOLERANCE;
    }
    return false;
}

bool LogVar::operator ==( lvFloat a ) const {
    LogVar other( a );
    return *this == other;
}

/* -------------------------------------------------------------- */

bool LogVar::operator !=( const LogVar& a ) const {
    return !( *this == a );
}

bool LogVar::operator !=( lvFloat a ) const {
    LogVar other( a );
    return *this != other;
}

/* -------------------------------------------------------------- */

LogVar LogVar::operator +( const LogVar& a ) const {
    return LogVar( LogVar::_logAdd( _logVal, a._logVal ) );
}

LogVar LogVar::operator +( lvFloat a ) const {
    return *this + LogVar( a );
}

/* -------------------------------------------------------------- */

LogVar& LogVar::operator +=( const LogVar& a ) {
    _logVal = LogVar::_logAdd( _logVal, a._logVal );
    return *this;
}

LogVar& LogVar::operator +=( lvFloat a ) {
    _logVal = LogVar::_logAdd( _logVal,LogVar( a )._logVal );
    return *this;
}

/* -------------------------------------------------------------- */

LogVar LogVar::operator *( const LogVar& a ) const {
    return LogVar( LogVar::_logMul( _logVal, a._logVal ) );
}

LogVar LogVar::operator *( lvFloat a ) const {
    return *this * LogVar( a );
}

/* -------------------------------------------------------------- */

LogVar& LogVar::operator *=( const LogVar& a ) {
    _logVal = LogVar::_logMul( _logVal, a._logVal );
    return *this;
}

LogVar& LogVar::operator *=( lvFloat a ) {
    _logVal = LogVar::_logMul( _logVal,LogVar( a )._logVal );
    return *this;
}

/* -------------------------------------------------------------- */

LogVar LogVar::operator /( const LogVar& a ) const {
    return LogVar( LogVar::_logDiv( _logVal, a._logVal ) );
}

LogVar LogVar::operator /( lvFloat a ) const {
    return *this / LogVar( a );
}

/* -------------------------------------------------------------- */

LogVar& LogVar::operator /=( const LogVar& a ) {
    _logVal = LogVar::_logDiv( _logVal, a._logVal );
    return *this;
}

LogVar& LogVar::operator /=( lvFloat a ) {
    _logVal = LogVar::_logDiv( _logVal,LogVar( a )._logVal );
    return *this;
}

/* -------------------------------------------------------------- */

lvFloat LogVar::toFloat() const {
    if( _logVal.i == ZERO ) {
        return 0.0;
    }
    return exp( _logVal.f );
}

/* -------------------------------------------------------------- */

LogVar operator +( lvFloat a, const LogVar& b ) {
    return b + a;
}

LogVar operator *( lvFloat a, const LogVar& b ) {
    return b * a;
}

LogVar operator /( lvFloat a, const LogVar& b ) {
    return b / a;
}
std::ostream& operator <<(std::ostream& os, const LogVar& lv) {
    os << lv.toFloat();
    return os;
}








