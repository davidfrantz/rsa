/**+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Named constant definitions
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**/


#ifndef CONSTANT_H
#define CONSTANT_H

#ifdef __cplusplus
extern "C" {
#endif

// string length
#define STRLEN 1024

// function return codes
#define SUCCESS 0
#define FAILURE 1
#define CANCEL 9

enum { NPOW_00 = 1,    NPOW_01 = 2,     NPOW_02 = 4,     NPOW_03 = 8, 
       NPOW_04 = 16,   NPOW_05 = 32,    NPOW_06 = 64,    NPOW_07 = 128,   
       NPOW_08 = 256,  NPOW_09 = 512,   NPOW_10 = 1024,  NPOW_11 = 2048, 
       NPOW_12 = 4096, NPOW_13 = 8192,  NPOW_14 = 16384, NPOW_15 = 32768, 
       NPOW_16 = 65536 };

#ifdef __cplusplus
}
#endif

#endif

