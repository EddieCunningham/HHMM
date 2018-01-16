#include "LogVar.h"
#include <cassert>
#include <cmath>

// https://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf
lvFloat log1pExp( lvFloat f ) {
    if( x < -37 ) {
        return exp( f );
    }
    else if( x <= 18 ) {
        return log1p( exp( f ) );
    }
    else if( x <= 33.3 ) {
        return f + exp( -f )
    }
    else {
        return f;
    }
}

lvFloat log1mExp( lvFloat f ) {
    if( f > -log( 2 ) ) {
        return log( -( exp( f ) - 1 ) )
    }
    return log1p( -( exp( f ) ) )
}

lvFloat logAdd( lvFloat logA, lvFloat logB ) {
    return logA + log1pExp( logB - logA );
}

lvFloat _logSub( lvFloat logA, lvFloat logB ) {
    return logA + log1mExp( logB - logA );
}

/* -------------------------------------------------------------- */

lv LogVar::logAdd( lv logX, lv logY ) {

    lvFloat ans;

    if( logX.i == ZERO ) {
        ans = logY;
    }
    if( logY.i == ZERO ) {
        ans = logX;
    }

    /* Don't do computations for adding what is practically 0 */
    if( logX.f - logY.f > LOG_UNDERFLOW ) {
        ans = logX;
    }
    if( logY.f - logX.f > LOG_UNDERFLOW ) {
        ans = logY;
    }

    if( logX.sgn == logY.sgn ) {
        // either computing ( x + y ) or -( x + y )
        if( logY.f > logX.f ) {
            ans.f = logAdd( logX, logY );
            ans.sgn = logX.sgn;
        }
        else {
            ans.f = logAdd( logY, logX );
            ans.sgn = logX.sgn;
        }
    }
    else {
        if( isClose( logX, logY ) ) {
            ans.i = ZERO;
        }
        else {
            // either computing ( a - b ) or -( a - b )
            if( logX.sgn == POSITIVE ) {
                // computing a - b
                if( logX.f > logY.f ) {
                    ans.f = _logSub( logX, logY );
                    ans.sgn = POSITIVE;
                }
                else {
                    ans.f = _logSub( logY, logX );
                    ans.sgn = NEGATIVE;
                }
            }
            else {
                // computing b - a
                if( logY.f > logX.f ) {
                    ans.f = _logSub( logY, logX );
                    ans.sgn = POSITIVE;
                }
                else {
                    ans.f = _logSub( logX, logY );
                    ans.sgn = NEGATIVE;
                }
            }
        }
    }

    return ans;
}

lv LogVar::logMul( lv logX, lv logY ) {

    lv ans;
    if( logX.i == ZERO || logY.i == ZERO ) {
        ans.i = ZERO;
    }
    else {
        ans.f = logX.f + logY.f;
        if( logX.sgn == logY.sgn ) {
            ans.sgn = POSITIVE;
        }
        else {
            ans.sgn = NEGATIVE;
        }
    }
    return ans;
}

lv LogVar::logDiv( lv logX, lv logY ) {

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
        if( logX.sgn == logY.sgn ) {
            ans.sgn = POSITIVE;
        }
        else {
            ans.sgn = NEGATIVE;
        }
    }
    return ans;
}

lv LogVar::parseFloat( lvFloat f, bool useSgn, bool sgn ) {

    /* Catch case where f == 0 */
    lv ans;
    if( FP_ZERO == fpclassify( f ) ) {
        ans.i = ZERO;
    }
    else {
        if( useSgn ) {
            if( f > 0 ) {
                ans.f = log( f );
            }
            else {
                ans.f = log( -f );
            }
            ans.sgn = sgn;
        }
        else {
            if( ans.f > 0 ) {
                ans.f = log( f );
                ans.sgn = POSITIVE;
            }
            else {
                ans.f = log( -f );
                ans.sgn = NEGATIVE;
            }
        }
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
    _logVal = LogVar::_parseFloat( a, false, true );
}

LogVar::LogVar( lvFloat a, bool alreadyLog ) {
    if( alreadyLog ) {
        if( FP_ZERO == fpclassify( a ) ) {
            _logVal.sgn = POSITIVE;
            _logVal.i = ZERO;
        }
        else if( a > 0 ) {
            _logVal.f = a;
            _logVal.sgn = POSITIVE;
            _logVal.i = NOTZERO;
        }
        else if( a < 0 ) {
            _logVal.f = -a;
            _logVal.sgn = NEGATIVE;
            _logVal.i = NOTZERO;
        }
    }
    else {
        /* Safely parse a float */
        _logVal = LogVar::_parseFloat( a, false, true );
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
    return LogVar( LogVar::logAdd( _logVal, a._logVal ) );
}

LogVar LogVar::operator +( lvFloat a ) const {
    return *this + LogVar( a );
}

/* -------------------------------------------------------------- */

LogVar& LogVar::operator +=( const LogVar& a ) {
    _logVal = LogVar::logAdd( _logVal, a._logVal );
    return *this;
}

LogVar& LogVar::operator +=( lvFloat a ) {
    _logVal = LogVar::logAdd( _logVal, LogVar( a )._logVal );
    return *this;
}

/* -------------------------------------------------------------- */

LogVar LogVar::operator -( LogVar a ) const {
    if( a.sgn == POSITIVE ) {
        a.sgn = NEGATIVE;
    }
    else {
        a.sgn = POSITIVE;
    }
    return LogVar( LogVar::logAdd( _logVal, a._logVal ) );
}

LogVar LogVar::operator -( lvFloat a ) const {
    LogVar ans = *this - LogVar( a );
    if( ans.sgn == POSITIVE ) {
        ans.sgn = NEGATIVE;
    }
    else {
        ans.sgn = POSITIVE;
    }
    return ans;
}

/* -------------------------------------------------------------- */

LogVar& LogVar::operator -=( LogVar a ) {
    if( a.sgn == POSITIVE ) {
        a.sgn = NEGATIVE;
    }
    else {
        a.sgn = POSITIVE;
    }
    _logVal = LogVar::logAdd( _logVal, a._logVal );
    return *this;
}

LogVar& LogVar::operator -=( lvFloat a ) {
    LogVar copy = LogVar( a );
    if( copy._logVal.sgn == POSITIVE ) {
        copy._logVal.sgn = NEGATIVE;
    }
    else {
        copy._logVal.sgn = POSITIVE;
    }
    _logVal = LogVar::logAdd( _logVal, copy );
    return *this;
}

/* -------------------------------------------------------------- */

LogVar LogVar::operator *( const LogVar& a ) const {
    return LogVar( LogVar::logMul( _logVal, a._logVal ) );
}

LogVar LogVar::operator *( lvFloat a ) const {
    return *this * LogVar( a );
}

/* -------------------------------------------------------------- */

LogVar& LogVar::operator *=( const LogVar& a ) {
    _logVal = LogVar::logMul( _logVal, a._logVal );
    return *this;
}

LogVar& LogVar::operator *=( lvFloat a ) {
    _logVal = LogVar::logMul( _logVal, LogVar( a )._logVal );
    return *this;
}

/* -------------------------------------------------------------- */

LogVar LogVar::operator /( const LogVar& a ) const {
    return LogVar( LogVar::logDiv( _logVal, a._logVal ) );
}

LogVar LogVar::operator /( lvFloat a ) const {
    return *this / LogVar( a );
}

/* -------------------------------------------------------------- */

LogVar& LogVar::operator /=( const LogVar& a ) {
    _logVal = LogVar::logDiv( _logVal, a._logVal );
    return *this;
}

LogVar& LogVar::operator /=( lvFloat a ) {
    _logVal = LogVar::logDiv( _logVal, LogVar( a )._logVal );
    return *this;
}

/* -------------------------------------------------------------- */

lvFloat LogVar::toFloat() const {
    if( _logVal.i == ZERO ) {
        return 0.0;
    }
    if( _logVal.sgn == POSITIVE ) {
        return exp( _logVal.f );
    }
    else {
        return -exp( _logVal.f );
    }
}

lvFloat LogVar::logValue() const {
    if( _logVal.i == ZERO ) {
        assert( 0 );
        return 0;
    }
    if( _logVal.sgn == POSITIVE ) {
        return _logVal.f;
    }
    else {
        return -_logVal.f;
    }
}

/* -------------------------------------------------------------- */

LogVar operator +( lvFloat a, const LogVar& b ) {
    return b + a;
}

LogVar operator -( lvFloat a, const LogVar& b ) {
    return b - a;
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








