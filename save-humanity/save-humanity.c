#include <stdio.h>
#include <string.h>

#define MAX_LEN 100000
#define K 1u
#define ALPHABET_SIZE 26

#define xstr(x) str(x)
#define str(x) #x

static unsigned int dk[K+1][ALPHABET_SIZE];

unsigned int umin( unsigned int x, unsigned int y )
{
  return (x < y) ? x : y;
}

unsigned int umax( unsigned int x, unsigned int y )
{
  return (x > y) ? x : y;
}
  
int read_num_cases( unsigned int *num_cases )
{
  return scanf( "%d", num_cases ) != 1;
}

int read_case( char *patient, char *virus )
{
  return scanf( "%s", patient ) != 1 || scanf( "%s", virus ) != 1;
}


int fill_dk( const char *P, unsigned int m )
{
  unsigned int ready[ ALPHABET_SIZE ];
  unsigned int a, i;

  for ( a = 0; a < ALPHABET_SIZE; ++a ) {
    ready[ a ] = m + 1;
  }

  for ( a = 0; a < ALPHABET_SIZE; ++a ) {
    for ( i = 0; i < K + 1; ++i ) {
#ifdef DEBUG
      if ( i >= K + 1 || a >= ALPHABET_SIZE ) {
        fprintf( stderr, "Error. Index out of bounds (1).\n" );
        return 1;
      }
#endif
      dk[i][a] = m;
    }
  }

  for ( i = m - 1; i >= 1; --i ) {
    unsigned int j, end;
    unsigned int pi = P[i - 1] - 'a';
    end = umax( i + 1, m - K );
    for ( j = ready[ pi ] - 1; j >= end; --j ) {
#ifdef DEBUG
      if ( j + K - m >= K + 1 || pi >= ALPHABET_SIZE ) {
        fprintf( stderr, "Error. Index out of bounds (2).\n" );
        return 1;
      }
#endif
      dk[j + K - m][pi] = j - i;
    }

    ready[pi] = umax( i + 1, m - K);
  }

#ifdef DEBUG
  printf( " i" );
  for ( i = 0; i < ALPHABET_SIZE; ++i ) {
    printf( "  %c", 'a' + i );
  }
  printf( "\n" );
  for ( i = 0; i < K + 1; ++i ) {
    int j;
    printf( "%2d", m - K + i );
    for ( j = 0; j < ALPHABET_SIZE; ++j ) {
      printf( " %2d", dk[i][j] );
    }
    printf( "\n" );
  }
#endif
  return 0;
}

int abm_search( int *out_matches, const char *T, const char *P, 
    unsigned int m )
{
  unsigned int n, j, o;
  n = strlen( T );
  o = 0;

  j = m;

  while ( j <= n ) {
    unsigned int h, i, neq, d;
    h = j; i = m; neq = 0;
    d = m - K;

    while ( i > 0 && neq <= K ) {
      unsigned int th, pi;
      th = T[ h-1 ] - 'a';

      if ( i >= m - K ) { 
        unsigned int dk_i_th; 
        dk_i_th = dk[i + K - m][ th ];
        d = umin( d, dk_i_th );
      }

      pi = P[i-1] - 'a';
      if ( th != pi ) {
        ++neq;
      }

      --i;
      --h;
    }

    if ( neq <= K ) {
      out_matches[o++] = j - m;
    }

    j += d;
  }

  out_matches[o] = -1;

  return 0;
}

int find_matches( int *out_matches, const char *text, const char *word )
{
  /* This function implements the Approximate Boyer-Moore k-mismatch string
   * matching described in:
   *
   * http://www.cs.hut.fi/~tarhio/papers/abm.pdf
   */

  int error = 0;
  unsigned int word_len;

  word_len = strlen( word );

  /* handle the special case where the word is one character long
   */
  if ( word_len <= K ) {
    unsigned int text_len, matches, i;
    text_len = strlen( text );
    matches = text_len - word_len + 1;
    for ( i = 0; i < matches; ++i ) {
      out_matches[i] = i;
    }
    out_matches[matches] = -1;
    return 0;
  }

  error = fill_dk( word, word_len );
  if ( error ) {
    fprintf( stderr, "Error making dk table.\n" );
    return error;
  }

  error = abm_search( out_matches, text, word, word_len );
  if ( error ) {
    fprintf( stderr, "Error making dk table.\n" );
    return error;
  }

  return 0;
}

int print_matches( const int *matches )
{
  int i;

  if ( matches[0] < 0 ) {
    printf( "\n" );
    return 0;
  }

  printf( "%d", matches[0] );

  for ( i = 1; matches[i] >= 0; ++i ) {
    printf( " %d", matches[i] );
  }
  printf( "\n" );

  return 0;
}

int main()
{
  int error;
  unsigned int num_cases, i;
  char patient[MAX_LEN + 1], virus[MAX_LEN + 1];
  int matches[MAX_LEN +1];

  error = read_num_cases( &num_cases );
  if ( error ) {
    fprintf( stderr, "Error reading number of cases.\n" );
    return error;
  }

#ifdef DEBUG
  printf( "Number of cases: %d\n", num_cases );
#endif

  for ( i = 0; i < num_cases; ++i ) {

    error = read_case( &patient[0], &virus[0] );
    if ( error ) {
      fprintf( stderr, "Error reading case %d\n", i );
      return error;
    }

#ifdef DEBUG
    printf( "Case %d:\n", i );
    printf( "  Patient: %s\n", &patient[0] );
    printf( "  Virus:   %s\n", &virus[0] );
#endif

    error = find_matches( &matches[0], &patient[0], &virus[0] );
    if ( error ) {
      fprintf( stderr, "Error finding matches in case %d\n", i );
      return error;
    }

    error = print_matches( &matches[0] );
    if ( error ) {
      fprintf( stderr, "Error printing matches in case %d\n", i );
      return error;
    }
  }

  return 0;
}

/* vim:se ts=2 sw=2 et: */

