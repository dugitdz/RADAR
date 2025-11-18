// serial_to_csv_sync_defs.c (com print no terminal)
// Lê a porta serial no Windows e grava dois CSVs sincronizados,
// mantendo os dados antigos e mostrando tudo no terminal.

#define _CRT_SECURE_NO_WARNINGS
#include <windows.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

// ======== CONFIGURAÇÃO POR DEFINE ========
#define PORT_NAME   "COM5"
#define BAUD_RATE   115200
// ========================================

#define READ_BUF_SZ  4096
#define LINE_BUF_SZ  1024

static bool set_comm_state(HANDLE h, DWORD baud){
    DCB dcb = {0};
    dcb.DCBlength = sizeof(DCB);
    if(!GetCommState(h,&dcb)) return false;
    dcb.BaudRate = baud;
    dcb.ByteSize = 8;
    dcb.Parity   = NOPARITY;
    dcb.StopBits = ONESTOPBIT;
    dcb.fOutxCtsFlow = FALSE;
    dcb.fOutxDsrFlow = FALSE;
    dcb.fDtrControl  = DTR_CONTROL_DISABLE;
    dcb.fOutX        = FALSE;
    dcb.fInX         = FALSE;
    dcb.fRtsControl  = RTS_CONTROL_DISABLE;
    return SetCommState(h,&dcb)!=0;
}

static bool set_timeouts(HANDLE h){
    COMMTIMEOUTS t={0};
    t.ReadIntervalTimeout         = 50;
    t.ReadTotalTimeoutConstant    = 50;
    t.ReadTotalTimeoutMultiplier  = 0;
    t.WriteTotalTimeoutConstant   = 50;
    t.WriteTotalTimeoutMultiplier = 0;
    return SetCommTimeouts(h,&t)!=0;
}

static bool extract_value_after_colon(const char* line, const char* keyword, float* out){
    const char* k = strstr(line, keyword);
    if(!k) return false;
    const char* colon = strchr(k, ':');
    if(!colon) return false;
    const char* p = colon + 1;
    while(*p==' '||*p=='\t') ++p;
    char* endp = NULL;
    double v = strtod(p, &endp);
    if(endp == p) return false;
    *out = (float)v;
    return true;
}

static void make_port_name(const char* in, char* out, size_t outsz){
    if (_strnicmp(in, "\\\\.\\", 4) == 0) {
        snprintf(out, outsz, "%s", in);
        return;
    }
    snprintf(out, outsz, "\\\\.\\%s", in);
}

int main(void){
    // stdout sem buffer para imprimir em tempo real
    setvbuf(stdout, NULL, _IONBF, 0);

    char portName[64];
    make_port_name(PORT_NAME, portName, sizeof(portName));

    HANDLE h = CreateFileA(
        portName,
        GENERIC_READ | GENERIC_WRITE,
        0, NULL, OPEN_EXISTING,
        FILE_ATTRIBUTE_NORMAL,
        NULL
    );
    if(h == INVALID_HANDLE_VALUE){
        fprintf(stderr, "Erro ao abrir %s\n", portName);
        return 1;
    }

    if(!SetupComm(h, READ_BUF_SZ, READ_BUF_SZ))
        fprintf(stderr, "Aviso: SetupComm falhou (seguindo assim mesmo).\n");
    if(!set_comm_state(h, (DWORD)BAUD_RATE)){
        fprintf(stderr, "Erro configurando DCB (baud/paridade/stopbits).\n");
        CloseHandle(h);
        return 2;
    }
    if(!set_timeouts(h)){
        fprintf(stderr, "Erro configurando timeouts.\n");
        CloseHandle(h);
        return 3;
    }

    FILE* fHeart = fopen("heart_rate.csv","a");
    FILE* fBreath= fopen("breath_rate.csv","a");
    if(!fHeart || !fBreath){
        fprintf(stderr, "Erro abrindo CSVs para escrita.\n");
        if(fHeart) fclose(fHeart);
        if(fBreath) fclose(fBreath);
        CloseHandle(h);
        return 4;
    }

    // Cabeçalho se arquivo vazio
    fseek(fHeart, 0, SEEK_END);
    long szHeart = ftell(fHeart);
    fseek(fBreath, 0, SEEK_END);
    long szBreath = ftell(fBreath);

    if (szHeart == 0) fprintf(fHeart, "dt_seconds,value\n");
    if (szBreath == 0) fprintf(fBreath,"dt_seconds,value\n");

    // Marca novo dataset com data/hora
    time_t now_t = time(NULL);
    struct tm* tm_info = localtime(&now_t);
    char timestamp[64];
    strftime(timestamp, sizeof(timestamp), "%Y-%m-%d %H:%M:%S", tm_info);
    fprintf(fHeart, "# --- NOVO DATASET %s ---\n", timestamp);
    fprintf(fBreath,"# --- NOVO DATASET %s ---\n", timestamp);
    fflush(fHeart); fflush(fBreath);

    setvbuf(fHeart, NULL, _IOLBF, 0);
    setvbuf(fBreath, NULL, _IOLBF, 0);

    LARGE_INTEGER freq; QueryPerformanceFrequency(&freq);
    LARGE_INTEGER t0;
    bool have_t0 = false;

    char lineBuf[LINE_BUF_SZ]; size_t lineLen=0;
    unsigned char rbuf[READ_BUF_SZ];

    printf("Lendo %s @ %d baud...\n", portName, (int)BAUD_RATE);
    printf("t0 = primeira 'Distance:' recebida (fallback: primeira linha válida).\n");
    printf("Arquivos: heart_rate.csv, breath_rate.csv\n");
    printf("===============================================\n");

    for(;;){
        DWORD nread=0;
        if(!ReadFile(h, rbuf, sizeof(rbuf), &nread, NULL)){
            fprintf(stderr,"Erro de leitura na serial.\n");
            break;
        }
        if(nread==0) continue;

        for(DWORD i=0;i<nread;++i){
            char c=(char)rbuf[i];
            if(c=='\r') continue;
            if(c=='\n'){
                lineBuf[lineLen]='\0';

                LARGE_INTEGER now; QueryPerformanceCounter(&now);

                float val=0.f;
                if(!have_t0 && extract_value_after_colon(lineBuf, "Distance", &val)){
                    t0 = now; have_t0 = true;
                }
                if(!have_t0){
                    float tmp=0.f;
                    if( extract_value_after_colon(lineBuf,"Heart Rate",&tmp) ||
                        extract_value_after_colon(lineBuf,"Breath Rate",&tmp) ||
                        extract_value_after_colon(lineBuf,"Distance",&tmp) ){
                        t0 = now; have_t0 = true;
                    }
                }

                // Heart Rate
                if(extract_value_after_colon(lineBuf, "Heart Rate", &val)){
                    if(val > 0.f && have_t0){
                        double dt = (double)(now.QuadPart - t0.QuadPart) / (double)freq.QuadPart;
                        fprintf(fHeart, "%.6f,%.2f\n", dt, val);
                        printf("[Heart] t=%.3fs -> %.2f bpm\n", dt, val);
                    }
                }
                // Breath Rate
                else if(extract_value_after_colon(lineBuf, "Breath Rate", &val)){
                    if(val > 0.f && have_t0){
                        double dt = (double)(now.QuadPart - t0.QuadPart) / (double)freq.QuadPart;
                        fprintf(fBreath, "%.6f,%.2f\n", dt, val);
                        printf("[Breath] t=%.3fs -> %.2f bpm\n", dt, val);
                    }
                }

                lineLen = 0;
            } else {
                if(lineLen + 1 < LINE_BUF_SZ) {
                    lineBuf[lineLen++] = c;
                } else {
                    lineLen = 0;
                }
            }
        }
    }

    fclose(fHeart);
    fclose(fBreath);
    CloseHandle(h);
    return 0;
}
