// 3_fourier.c
// Serie de Fourier con validaciones robustas

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <errno.h>

// ============================================================================
// PARAMETROS CONFIGURABLES
// ============================================================================
#define FUNCION_ORIGINAL(x) ((x) < L ? (x) : (2*L - (x)))
#define L                   M_PI
#define N_TERMINOS          10
#define PUNTOS_GRAFICO      500
#define GRAFICO_INICIO      0.0
#define GRAFICO_FIN         2*M_PI
#define NOMBRE_GRAFICO      "fourier_grafico.png"
#define ANCHO_GRAFICO       800
#define ALTO_GRAFICO        600
// ============================================================================

// ============================================================================
// FUNCIONES DE VALIDACION
// ============================================================================
int es_numerico_valido(double valor) {
    return !(isnan(valor) || isinf(valor) || fabs(valor) > 1e100);
}

void verificar_nan_inf(const char *nombre, double valor, int linea) {
    if (isnan(valor)) {
        printf("ERROR en linea %d: %s = NaN\n", linea, nombre);
        exit(EXIT_FAILURE);
    }
    if (isinf(valor)) {
        printf("ERROR en linea %d: %s = Infinito\n", linea, nombre);
        exit(EXIT_FAILURE);
    }
}

#define VALIDAR(variable) verificar_nan_inf(#variable, variable, __LINE__)

FILE* abrir_archivo(const char *nombre, const char *modo) {
    FILE *archivo = fopen(nombre, modo);
    if (archivo == NULL) {
        printf("ERROR: No se pudo abrir '%s'\n", nombre);
        exit(EXIT_FAILURE);
    }
    return archivo;
}

// ============================================================================
// VALIDACION DE PARAMETROS
// ============================================================================
void validar_parametros() {
    printf("Validando parametros...\n");
    
    if (L <= 0) {
        printf("ERROR: L debe ser positivo (L = %f)\n", L);
        exit(EXIT_FAILURE);
    }
    
    if (N_TERMINOS <= 0) {
        printf("ERROR: N_TERMINOS debe ser positivo (%d)\n", N_TERMINOS);
        exit(EXIT_FAILURE);
    }
    
    if (PUNTOS_GRAFICO < 10) {
        printf("ERROR: PUNTOS_GRAFICO debe ser >= 10 (%d)\n", PUNTOS_GRAFICO);
        exit(EXIT_FAILURE);
    }
    
    if (GRAFICO_FIN <= GRAFICO_INICIO) {
        printf("ERROR: GRAFICO_FIN debe ser > GRAFICO_INICIO\n");
        exit(EXIT_FAILURE);
    }
    
    // Validar funcion en algunos puntos
    for (int i = 0; i < 5; i++) {
        double x = GRAFICO_INICIO + i * (GRAFICO_FIN - GRAFICO_INICIO) / 4;
        double fx = FUNCION_ORIGINAL(x);
        VALIDAR(fx);
        
        if (!es_numerico_valido(fx)) {
            printf("ERROR: Funcion invalida en x = %f\n", x);
            exit(EXIT_FAILURE);
        }
    }
    
    printf("Parametros validados correctamente\n\n");
}

// ============================================================================
// FUNCIONES PRINCIPALES
// ============================================================================
int main() {
    validar_parametros();
    
    printf("==============================================================\n");
    printf("                    SERIE DE FOURIER                          \n");
    printf("==============================================================\n\n");
    
    printf("CONFIGURACION:\n");
    printf("  Funcion:          Triangular en [0, 2pi]\n");
    printf("  Periodo:          L = pi\n");
    printf("  Terminos:         %d\n", N_TERMINOS);
    printf("  Puntos grafico:   %d\n\n", PUNTOS_GRAFICO);
    
    // ============================================================================
    // CALCULAR COEFICIENTES CON VALIDACION
    // ============================================================================
    printf("CALCULANDO COEFICIENTES...\n");
    printf("-------------------------------------------------------------\n");
    
    double a0 = 0.0;
    double an[N_TERMINOS + 1], bn[N_TERMINOS + 1];
    
    int puntos_int = 1000;
    double dx_int = 2*L / puntos_int;
    VALIDAR(dx_int);
    
    if (dx_int <= 0) {
        printf("ERROR: dx_int no valido: %f\n", dx_int);
        return EXIT_FAILURE;
    }
    
    // Calcular a0
    for (int i = 0; i < puntos_int; i++) {
        double x = i * dx_int;
        double f = FUNCION_ORIGINAL(x);
        VALIDAR(f);
        
        if (!es_numerico_valido(f)) {
            printf("ERROR: Funcion invalida en x = %f, f(x) = %f\n", x, f);
            return EXIT_FAILURE;
        }
        
        a0 += f;
        VALIDAR(a0);
    }
    a0 *= dx_int / L;
    VALIDAR(a0);
    
    printf("  Coeficiente a0 = %.6f\n", a0);
    
    // Calcular coeficientes an y bn
    for (int n = 1; n <= N_TERMINOS; n++) {
        double suma_an = 0.0, suma_bn = 0.0;
        
        for (int i = 0; i < puntos_int; i++) {
            double x = i * dx_int;
            double f = FUNCION_ORIGINAL(x);
            VALIDAR(f);
            
            double cos_val = cos(n * M_PI * x / L);
            double sin_val = sin(n * M_PI * x / L);
            VALIDAR(cos_val); VALIDAR(sin_val);
            
            suma_an += f * cos_val;
            suma_bn += f * sin_val;
            
            VALIDAR(suma_an); VALIDAR(suma_bn);
            
            // Detectar overflow
            if (fabs(suma_an) > 1e50 || fabs(suma_bn) > 1e50) {
                printf("ERROR: Overflow en calculo de coeficientes n=%d\n", n);
                return EXIT_FAILURE;
            }
        }
        
        an[n] = suma_an * dx_int / L;
        bn[n] = suma_bn * dx_int / L;
        
        VALIDAR(an[n]); VALIDAR(bn[n]);
        
        if (n <= 5) {
            printf("  a%d = %9.6f, b%d = %9.6f\n", n, an[n], n, bn[n]);
        }
    }
    
    // ============================================================================
    // GENERAR DATOS CON VALIDACION
    // ============================================================================
    printf("\nGENERANDO DATOS...\n");
    printf("-------------------------------------------------------------\n");
    
    FILE *orig = abrir_archivo("fourier_original.dat", "w");
    FILE *serie = abrir_archivo("fourier_serie.dat", "w");
    FILE *script_gp = abrir_archivo("fourier_plot.gp", "w");
    
    fprintf(orig, "# Funcion original\n");
    fprintf(serie, "# Aproximacion de Fourier\n");
    
    double dx = (GRAFICO_FIN - GRAFICO_INICIO) / PUNTOS_GRAFICO;
    VALIDAR(dx);
    
    int errores_puntos = 0;
    
    for (int i = 0; i <= PUNTOS_GRAFICO; i++) {
        double x = GRAFICO_INICIO + i * dx;
        VALIDAR(x);
        
        // Funcion original
        double f_orig = FUNCION_ORIGINAL(x);
        VALIDAR(f_orig);
        
        // Serie de Fourier
        double f_serie = a0 / 2;
        VALIDAR(f_serie);
        
        for (int n = 1; n <= N_TERMINOS; n++) {
            double termino = an[n] * cos(n * M_PI * x / L) + bn[n] * sin(n * M_PI * x / L);
            VALIDAR(termino);
            
            f_serie += termino;
            VALIDAR(f_serie);
            
            // Detectar divergencia
            if (!es_numerico_valido(f_serie)) {
                printf("ADVERTENCIA: Serie divergente en x=%.3f, n=%d\n", x, n);
                errores_puntos++;
                f_serie = 0; // Reset para continuar
                break;
            }
        }
        
        // Guardar si ambos valores son validos
        if (es_numerico_valido(f_orig) && es_numerico_valido(f_serie)) {
            fprintf(orig, "%.6f %.6f\n", x, f_orig);
            fprintf(serie, "%.6f %.6f\n", x, f_serie);
        } else {
            errores_puntos++;
        }
        
        // Progreso
        if (i % (PUNTOS_GRAFICO/10) == 0) {
            printf("  %3d%%: x=%.3f, f(x)=%.3f, Fourier=%.3f\n", 
                   (i*100)/PUNTOS_GRAFICO, x, f_orig, f_serie);
        }
    }
    
    fclose(orig);
    fclose(serie);
    
    if (errores_puntos > 0) {
        printf("ADVERTENCIA: %d puntos tuvieron problemas numericos\n", errores_puntos);
    }
    
    // ============================================================================
    // CREAR SCRIPT GNUPLOT
    // ============================================================================
    fprintf(script_gp, "# Script para serie de Fourier\n");
    fprintf(script_gp, "set terminal pngcairo size %d,%d enhanced font 'Arial,10'\n", 
            ANCHO_GRAFICO, ALTO_GRAFICO);
    fprintf(script_gp, "set output '%s'\n", NOMBRE_GRAFICO);
    fprintf(script_gp, "set title 'Serie de Fourier (N = %d terminos)'\n", N_TERMINOS);
    fprintf(script_gp, "set xlabel 'x'\n");
    fprintf(script_gp, "set ylabel 'f(x)'\n");
    fprintf(script_gp, "set grid\n");
    fprintf(script_gp, "set key top left box\n");
    fprintf(script_gp, "set xrange [%f:%f]\n", GRAFICO_INICIO, GRAFICO_FIN);
    fprintf(script_gp, "set yrange [-0.5:4.5]\n\n");
    
    fprintf(script_gp, "plot 'fourier_original.dat' w l lw 3 lc rgb '#0066CC' title 'Funcion original', \\\n");
    fprintf(script_gp, "     'fourier_serie.dat' w l lw 2 lc rgb '#FF3333' dt 2 title 'Aproximacion Fourier'\n");
    
    fclose(script_gp);
    
    // ============================================================================
    // EJECUTAR GNUPLOT
    // ============================================================================
    printf("\nGENERANDO GRAFICO...\n");
    printf("-------------------------------------------------------------\n");
    
    int resultado = system("gnuplot fourier_plot.gp 2>&1");
    
    if (resultado != 0) {
        printf("ADVERTENCIA: Problema al generar grafico\n");
    } else {
        printf("Grafico generado: %s\n", NOMBRE_GRAFICO);
    }
    
    // ============================================================================
    // ANALISIS DE ERROR
    // ============================================================================
    printf("\nANALISIS DE ERROR:\n");
    printf("-------------------------------------------------------------\n");
    
    double error_cuadratico = 0.0;
    double error_maximo = 0.0;
    int puntos_error = 100;
    int puntos_validos = 0;
    
    for (int i = 0; i <= puntos_error; i++) {
        double x = GRAFICO_INICIO + i * (GRAFICO_FIN - GRAFICO_INICIO) / puntos_error;
        double f_orig = FUNCION_ORIGINAL(x);
        double f_serie = a0 / 2;
        
        for (int n = 1; n <= N_TERMINOS; n++) {
            f_serie += an[n] * cos(n * M_PI * x / L) + bn[n] * sin(n * M_PI * x / L);
        }
        
        if (es_numerico_valido(f_orig) && es_numerico_valido(f_serie)) {
            double error = fabs(f_orig - f_serie);
            error_cuadratico += error * error;
            if (error > error_maximo) error_maximo = error;
            puntos_validos++;
        }
    }
    
    if (puntos_validos > 0) {
        error_cuadratico = sqrt(error_cuadratico / puntos_validos);
        printf("  Error cuadratico medio: %.6f\n", error_cuadratico);
        printf("  Error maximo:           %.6f\n", error_maximo);
        printf("  Puntos analizados:      %d/%d\n", puntos_validos, puntos_error + 1);
    } else {
        printf("ERROR: No se pudieron calcular errores\n");
    }
    
    // ============================================================================
    // RESUMEN
    // ============================================================================
    printf("\nRESUMEN:\n");
    printf("-------------------------------------------------------------\n");
    printf("  Terminos calculados:  %d\n", N_TERMINOS);
    printf("  Puntos generados:     %d\n", PUNTOS_GRAFICO + 1);
    printf("  Errores encontrados:  %d\n", errores_puntos);
    printf("  Grafico:              %s\n", 
           (resultado == 0) ? "GENERADO" : "NO GENERADO");
    
    printf("\n==============================================================\n");
    printf("                      EJECUCION COMPLETADA                     \n");
    printf("==============================================================\n");
    
    return EXIT_SUCCESS;
}