// 4_derivadas.c
// Calculo de 8 derivadas numericas con validaciones robustas

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <errno.h>

// ============================================================================
// PARAMETROS CONFIGURABLES
// ============================================================================
#define FUNCION_X(x)        (sin(x) + (x)*(x))
#define FUNCION_XY(x,y)     ((x)*(x)*sin(y) + exp((x)*(y)))
#define PUNTO_X0            1.0
#define PUNTO_Y0            0.5
#define PASO_H              0.0001
#define GRAFICO_INICIO      (PUNTO_X0 - 2.0)
#define GRAFICO_FIN         (PUNTO_X0 + 2.0)
#define GRAFICO_PUNTOS      100
#define NOMBRE_GRAFICO      "derivadas_grafico.png"
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
        printf(" ERROR [Linea %d]: %s = NaN (Operacion invalida)\n", linea, nombre);
        printf("   Revise funciones matematicas (division por cero, raiz negativa, etc.)\n");
        exit(EXIT_FAILURE);
    }
    if (isinf(valor)) {
        printf(" ERROR [Linea %d]: %s = Infinito (Overflow)\n", linea, nombre);
        printf("   Valor demasiado grande, reduzca paso h o punto de evaluacion\n");
        exit(EXIT_FAILURE);
    }
    if (!es_numerico_valido(valor)) {
        printf(" ERROR [Linea %d]: %s = Valor numerico invalido\n", linea, nombre);
        exit(EXIT_FAILURE);
    }
}

#define VALIDAR(variable) verificar_nan_inf(#variable, variable, __LINE__)

void validar_parametro_h(double h) {
    if (h <= 0) {
        printf(" ERROR: Paso h debe ser positivo (h = %.6f)\n", h);
        exit(EXIT_FAILURE);
    }
    if (h > 1.0) {
        printf(" ADVERTENCIA: Paso h muy grande (h = %.6f)\n", h);
        printf("   Las derivadas pueden ser imprecisas\n");
    }
    if (h < 1e-10) {
        printf(" ADVERTENCIA: Paso h muy pequeño (h = %.6f)\n", h);
        printf("   Posible error de cancelacion numerica\n");
    }
}

FILE* abrir_archivo(const char *nombre, const char *modo) {
    FILE *archivo = fopen(nombre, modo);
    if (archivo == NULL) {
        printf(" ERROR: No se pudo abrir archivo '%s' (modo: %s)\n", nombre, modo);
        printf("   Codigo de error: %d\n", errno);
        printf("   Verifique permisos y espacio en disco\n");
        exit(EXIT_FAILURE);
    }
    return archivo;
}

void validar_funciones_punto(double x, double y) {
    double fx = FUNCION_X(x);
    double fxy = FUNCION_XY(x, y);
    
    if (!es_numerico_valido(fx)) {
        printf(" ERROR: f(x) invalida en x = %.6f\n", x);
        printf("   f(%.6f) = %.6f\n", x, fx);
        exit(EXIT_FAILURE);
    }
    
    if (!es_numerico_valido(fxy)) {
        printf(" ERROR: f(x,y) invalida en (%.6f, %.6f)\n", x, y);
        printf("   f(%.6f, %.6f) = %.6f\n", x, y, fxy);
        exit(EXIT_FAILURE);
    }
}

// ============================================================================
// FUNCIONES DE CALCULO CON VALIDACION
// ============================================================================
double calcular_derivada_primera(double x0, double h) {
    double f_plus = FUNCION_X(x0 + h);
    double f_minus = FUNCION_X(x0 - h);
    
    VALIDAR(f_plus);
    VALIDAR(f_minus);
    
    // Verificar diferencia significativa
    if (fabs(f_plus - f_minus) < 1e-15) {
        printf(" ADVERTENCIA: Diferencia muy pequeña en derivada primera\n");
        printf("   f(x+h) - f(x-h) = %.2e\n", f_plus - f_minus);
    }
    
    double derivada = (f_plus - f_minus) / (2*h);
    VALIDAR(derivada);
    
    return derivada;
}

double calcular_derivada_segunda(double x0, double h) {
    double f_plus = FUNCION_X(x0 + h);
    double f_center = FUNCION_X(x0);
    double f_minus = FUNCION_X(x0 - h);
    
    VALIDAR(f_plus);
    VALIDAR(f_center);
    VALIDAR(f_minus);
    
    // Verificar valores extremos
    double max_val = fmax(fabs(f_plus), fmax(fabs(f_center), fabs(f_minus)));
    if (max_val > 1e50) {
        printf(" ADVERTENCIA: Valores muy grandes en derivada segunda\n");
    }
    
    double derivada = (f_plus - 2*f_center + f_minus) / (h*h);
    VALIDAR(derivada);
    
    return derivada;
}

double calcular_derivada_parcial_x(double x0, double y0, double h) {
    double f_plus = FUNCION_XY(x0 + h, y0);
    double f_minus = FUNCION_XY(x0 - h, y0);
    
    VALIDAR(f_plus);
    VALIDAR(f_minus);
    
    double derivada = (f_plus - f_minus) / (2*h);
    VALIDAR(derivada);
    
    return derivada;
}

double calcular_derivada_parcial_y(double x0, double y0, double h) {
    double f_plus = FUNCION_XY(x0, y0 + h);
    double f_minus = FUNCION_XY(x0, y0 - h);
    
    VALIDAR(f_plus);
    VALIDAR(f_minus);
    
    double derivada = (f_plus - f_minus) / (2*h);
    VALIDAR(derivada);
    
    return derivada;
}

double calcular_derivada_mixta(double x0, double y0, double h) {
    double f_pp = FUNCION_XY(x0 + h, y0 + h);
    double f_pm = FUNCION_XY(x0 + h, y0 - h);
    double f_mp = FUNCION_XY(x0 - h, y0 + h);
    double f_mm = FUNCION_XY(x0 - h, y0 - h);
    
    VALIDAR(f_pp); VALIDAR(f_pm);
    VALIDAR(f_mp); VALIDAR(f_mm);
    
    // Verificar simetria (para funciones suaves)
    if (fabs(f_pp - f_mm) > 1e-6 * fmax(fabs(f_pp), fabs(f_mm))) {
        printf(" ADVERTENCIA: Asimetria en derivada mixta\n");
        printf("   f(x+h,y+h) - f(x-h,y-h) = %.2e\n", f_pp - f_mm);
    }
    
    double derivada = (f_pp - f_pm - f_mp + f_mm) / (4*h*h);
    VALIDAR(derivada);
    
    return derivada;
}

// ============================================================================
// PROGRAMA PRINCIPAL
// ============================================================================
int main() {
    // ============================================================================
    // VALIDACION INICIAL
    // ============================================================================
    printf(" VALIDANDO PARAMETROS INICIALES...\n");
    printf("-------------------------------------------------------------\n");
    
    validar_parametro_h(PASO_H);
    validar_funciones_punto(PUNTO_X0, PUNTO_Y0);
    
    double h = PASO_H;
    double x0 = PUNTO_X0, y0 = PUNTO_Y0;
    
    // Validar entorno de calculo
    if (!es_numerico_valido(x0) || !es_numerico_valido(y0)) {
        printf(" ERROR: Puntos de evaluacion invalidos\n");
        printf("   x0 = %.6f, y0 = %.6f\n", x0, y0);
        return EXIT_FAILURE;
    }
    
    printf(" Punto de evaluacion valido: (%.6f, %.6f)\n", x0, y0);
    printf(" Paso h valido: %.6f\n", h);
    printf(" Funciones validas en el punto\n\n");
    
    // ============================================================================
    // ENCABEZADO
    // ============================================================================
    printf(" 8 DERIVADAS NUMERICAS CON VALIDACION \n\n");
    
    printf("FUNCIONES:\n");
    printf("  f(x)   = sin(x) + x²\n");
    printf("  f(x,y) = x²·sin(y) + e^(x·y)\n");
    printf("PUNTO:   (x0, y0) = (%.1f, %.1f)\n", x0, y0);
    printf("PASO:    h = %.4f\n\n", h);
    
    // ============================================================================
    // CALCULO DE DERIVADAS CON VALIDACION
    // ============================================================================
    printf(" CALCULANDO DERIVADAS...\n");
    printf("-------------------------------------------------------------\n");
    
    printf("+----+--------------------------------------+-----------------+----------+\n");
    printf("| #  | Derivada                            | Valor Numerico  | Estado   |\n");
    printf("+----+--------------------------------------+-----------------+----------+\n");
    
    // Derivada a) D[f(x), x]
    double da = calcular_derivada_primera(x0, h);
    printf("| a) | D[f(x), x]                          | %14.6f | - VALIDO |\n", da);
    
    // Derivada b) D[f(x), {x, 2}]
    double db = calcular_derivada_segunda(x0, h);
    printf("| b) | D[f(x), {x, 2}]                     | %14.6f | - VALIDO |\n", db);
    
    // Derivada c) D[f(x,y), x]
    double dc = calcular_derivada_parcial_x(x0, y0, h);
    printf("| c) | D[f(x,y), x]                        | %14.6f | - VALIDO |\n", dc);
    
    // Derivada d) D[f(x,y), y]
    double dd = calcular_derivada_parcial_y(x0, y0, h);
    printf("| d) | D[f(x,y), y]                        | %14.6f | - VALIDO |\n", dd);
    
    // Derivada e) D[f(x,y), {x, 2}]
    double de = (FUNCION_XY(x0 + h, y0) - 2*FUNCION_XY(x0, y0) + FUNCION_XY(x0 - h, y0)) / (h*h);
    VALIDAR(de);
    printf("| e) | D[f(x,y), {x, 2}]                   | %14.6f | - VALIDO |\n", de);
    
    // Derivada f) D[f(x,y), {y, 2}]
    double df = (FUNCION_XY(x0, y0 + h) - 2*FUNCION_XY(x0, y0) + FUNCION_XY(x0, y0 - h)) / (h*h);
    VALIDAR(df);
    printf("| f) | D[f(x,y), {y, 2}]                   | %14.6f | - VALIDO |\n", df);
    
    // Derivada g) D[f(x,y), {x, y}]
    double dg = calcular_derivada_mixta(x0, y0, h);
    printf("| g) | D[f(x,y), {x, y}]                   | %14.6f | - VALIDO |\n", dg);
    
    // Derivada h) D[f(x,y), {x, 3}]
    double f_2h = FUNCION_XY(x0 + 2*h, y0);
    double f_h = FUNCION_XY(x0 + h, y0);
    double f_mh = FUNCION_XY(x0 - h, y0);
    double f_m2h = FUNCION_XY(x0 - 2*h, y0);
    
    VALIDAR(f_2h); VALIDAR(f_h);
    VALIDAR(f_mh); VALIDAR(f_m2h);
    
    double dh = (f_2h - 2*f_h + 2*f_mh - f_m2h) / (2*h*h*h);
    VALIDAR(dh);
    printf("| h) | D[f(x,y), {x, 3}]                   | %14.6f | - VALIDO |\n", dh);
    
    printf("+----+--------------------------------------+-----------------+----------+\n");
    
    // ============================================================================
    // GENERAR DATOS PARA GRAFICAS
    // ============================================================================
    printf("\n GENERANDO DATOS PARA GRAFICAS...\n");
    printf("-------------------------------------------------------------\n");
    
    FILE *datos = abrir_archivo("derivadas.dat", "w");
    fprintf(datos, "# x df/dx d^2f/dx^2\n");
    
    double dx_graf = (GRAFICO_FIN - GRAFICO_INICIO) / GRAFICO_PUNTOS;
    int puntos_validos = 0, puntos_invalidos = 0;
    
    for (int i = 0; i <= GRAFICO_PUNTOS; i++) {
        double x = GRAFICO_INICIO + i * dx_graf;
        
        // Calcular derivadas con validacion
        double d1, d2;
        int valido = 1;
        
        try_calc:
        d1 = (FUNCION_X(x + h) - FUNCION_X(x - h)) / (2*h);
        d2 = (FUNCION_X(x + h) - 2*FUNCION_X(x) + FUNCION_X(x - h)) / (h*h);
        
        if (!es_numerico_valido(d1) || !es_numerico_valido(d2)) {
            if (h > 1e-10) {
                // Intentar con h mas grande
                h *= 2;
                printf(" Ajustando h a %.2e para x = %.3f\n", h, x);
                goto try_calc;
            } else {
                puntos_invalidos++;
                valido = 0;
            }
        }
        
        if (valido) {
            fprintf(datos, "%.6f %.6f %.6f\n", x, d1, d2);
            puntos_validos++;
        }
        
        // Progreso
        if (i % (GRAFICO_PUNTOS/10) == 0 && valido) {
            printf("  %3d%%: x=%.3f, f'(x)=%.3f, f''(x)=%.3f\n", 
                   (i*100)/GRAFICO_PUNTOS, x, d1, d2);
        }
    }
    
    fclose(datos);
    
    if (puntos_invalidos > 0) {
        printf(" ADVERTENCIA: %d puntos no pudieron calcularse\n", puntos_invalidos);
        printf("   Se generaron %d puntos validos de %d\n", puntos_validos, GRAFICO_PUNTOS + 1);
    }
    
    // ============================================================================
    // CREAR SCRIPT GNUPLOT
    // ============================================================================
    FILE *script = abrir_archivo("derivadas_plot.gp", "w");
    
    fprintf(script, "# Script para derivadas numericas\n");
    fprintf(script, "set terminal pngcairo size %d,%d enhanced font 'Arial,10'\n", 
            ANCHO_GRAFICO, ALTO_GRAFICO);
    fprintf(script, "set output '%s'\n", NOMBRE_GRAFICO);
    fprintf(script, "set title 'Derivadas de f(x) = sin(x) + x²'\n");
    fprintf(script, "set xlabel 'x'\n");
    fprintf(script, "set ylabel 'Valor de derivada'\n");
    fprintf(script, "set grid\n");
    fprintf(script, "set key top left box\n");
    fprintf(script, "set xrange [%f:%f]\n", GRAFICO_INICIO, GRAFICO_FIN);
    
    fprintf(script, "plot 'derivadas.dat' u 1:2 w l lw 2 lc rgb '#0066CC' title 'Primera derivada f''(x)', \\\n");
    fprintf(script, "     '' u 1:3 w l lw 2 lc rgb '#FF3333' dt 2 title 'Segunda derivada f''''(x)', \\\n");
    fprintf(script, "     %f, %f w p pt 7 ps 2 lc rgb '#00AA00' title 'Punto (x0, %.3f)'\n", 
            x0, da, da);
    
    fclose(script);
    
    // ============================================================================
    // EJECUTAR GNUPLOT
    // ============================================================================
    printf("\n GENERANDO GRAFICO...\n");
    printf("-------------------------------------------------------------\n");
    
    int resultado = system("gnuplot derivadas_plot.gp 2>&1");
    
    if (resultado != 0) {
        printf(" ADVERTENCIA: Gnuplot reporto problemas\n");
        printf("   Comando: gnuplot derivadas_plot.gp\n");
        printf("   Verifique que Gnuplot este instalado correctamente\n");
    } else {
        printf(" Grafico generado: %s\n", NOMBRE_GRAFICO);
    }
    
    // ============================================================================
    // ANALISIS DE ERROR (comparacion analitica)
    // ============================================================================
    printf("\n ANALISIS DE ERROR (comparacion con valores analiticos):\n");
    printf("-------------------------------------------------------------\n");
    
    // Valores analiticos exactos en x0 = 1.0
    double da_analitica = cos(x0) + 2*x0;        // cos(1) + 2
    double db_analitica = -sin(x0) + 2;          // -sin(1) + 2
    
    double error_abs_a = fabs(da - da_analitica);
    double error_rel_a = 100 * error_abs_a / fabs(da_analitica);
    double error_abs_b = fabs(db - db_analitica);
    double error_rel_b = 100 * error_abs_b / fabs(db_analitica);
    
    printf("  Primera derivada (f'(x)):\n");
    printf("    Valor numerico:   %.8f\n", da);
    printf("    Valor analitico:  %.8f\n", da_analitica);
    printf("    Error absoluto:   %.2e\n", error_abs_a);
    printf("    Error relativo:   %.2e%%\n", error_rel_a);
    
    printf("\n  Segunda derivada (f''(x)):\n");
    printf("    Valor numerico:   %.8f\n", db);
    printf("    Valor analitico:  %.8f\n", db_analitica);
    printf("    Error absoluto:   %.2e\n", error_abs_b);
    printf("    Error relativo:   %.2e%%\n", error_rel_b);
    
    // Evaluar calidad de aproximacion
    printf("\n  EVALUACION DE LA APROXIMACION:\n");
    if (error_rel_a < 0.1 && error_rel_b < 0.1) {
        printf("    - Excelente precision (error < 0.1%%)\n");
    } else if (error_rel_a < 1.0 && error_rel_b < 1.0) {
        printf("    - Buena precision (error < 1%%)\n");
    } else if (error_rel_a < 5.0 && error_rel_b < 5.0) {
        printf("    - Precision aceptable (error < 5%%)\n");
    } else {
        printf("    - Precision pobre, considere ajustar h\n");
    }
    
    // ============================================================================
    // RESUMEN FINAL
    // ============================================================================
    printf("\n RESUMEN DE EJECUCION:\n");
    printf("-------------------------------------------------------------\n");
    printf("  Derivadas calculadas:  8/8 exitosas\n");
    printf("  Puntos para grafico:   %d/%d validos\n", puntos_validos, GRAFICO_PUNTOS + 1);
    printf("  Error maximo:          %.2e%%\n", fmax(error_rel_a, error_rel_b));
    printf("  Grafico generado:      %s\n", (resultado == 0) ? "SI" : "NO");
    printf("  Archivos creados:      derivadas.dat, derivadas_plot.gp\n");
    
    printf("\n EJECUCION COMPLETADA\n");
    
    return EXIT_SUCCESS;
}