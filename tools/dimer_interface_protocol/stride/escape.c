#include <stdarg.h>
#include <stdio.h>

int escape(int RetVal, char *format, ... ) {
va_list ptr;

va_start(ptr,format);
vfprintf(stderr,format,ptr);
va_end(ptr);
return(RetVal);
}

