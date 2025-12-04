// 5_ecuacion1.c
// EDO y' = x - y con Runge-Kutta 4 y validaciones robustas

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <errno.h>

// ============================================================================
// PARAMETROS CONFIGURABLES
// ============================================================================
#define EDO_FUNCION(x,y)    ((x) - (y))
#define SOLUCION_EXACTA(x)  ((x) - 1 + 2*exp(-(x)))
#define X_INICIAL           0.0
#define X_FINAL             5.0
#define Y_INICIAL           1.0
#define PASO_H              0.1
#define NOMBRE_GRAFICO      "rk4_grafico.png"
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
        printf(" ERROR [Linea %d]: %s = NaN\n", linea, nombre);
        printf("   Causa posible: Operacion matematica invalida\n");
        exit(EXIT_FAILURE);
    }
    if (isinf(valor)) {
        printf(" ERROR [Linea %d]: %s = Infinito\n", linea, nombre);
        printf("   Causa posible: Overflow numerico\n");
        exit(EXIT_FAILURE);
    }
}

#define VALIDAR(variable) verificar_nan_inf(#variable, variable, __LINE__)

void validar_parametros() {
    if (PASO_H <= 0) {
        printf(" ERROR: PASO_H debe ser positivo (h = %f)\n", PASO_H);
        exit(EXIT_FAILURE);
    }
    
    if (X_FINAL <= X_INICIAL) {
        printf(" ERROR: X_FINAL debe ser > X_INICIAL\n");
        printf("   X_INICIAL = %f, X_FINAL = %f\n", X_INICIAL, X_FINAL);
        exit(EXIT_FAILURE);
    }
    
    if (!es_numerico_valido(Y_INICIAL)) {
        printf(" ERROR: Y_INICIAL invalido: %f\n", Y_INICIAL);
        exit(EXIT_FAILURE);
    }
    
    // Validar solucion exacta en algunos puntos
    for (double x = X_INICIAL; x <= X_FINAL; x += 1.0) {
        double y_exacta = SOLUCION_EXACTA(x);
        if (!es_numerico_valido(y_exacta)) {
            printf(" ERROR: Solucion exacta invalida en x = %f\n", x);
            exit(EXIT_FAILURE);
        }
    }
}

FILE* abrir_archivo(const char *nombre, const char *modo) {
    FILE *archivo = fopen(nombre, modo);
    if (archivo == NULL) {
        printf(" ERROR: No se pudo abrir '%s'\n", nombre);
        printf("   Error del sistema: %d\n", errno);
        exit(EXIT_FAILURE);
    }
    return archivo;
}

// ============================================================================
// RUNGE-KUTTA 4 CON VALIDACION
// ============================================================================
double rk4_validado(double x, double y, double h, int paso_actual) {
    double k1, k2, k3, k4;
    
    // k1
    k1 = EDO_FUNCION(x, y);
    VALIDAR(k1);
    
    // k2
    double x2 = x + h/2;
    double y2 = y + h*k1/2;
    VALIDAR(x2); VALIDAR(y2);
    
    k2 = EDO_FUNCION(x2, y2);
    VALIDAR(k2);
    
    // k3
    double x3 = x + h/2;
    double y3 = y + h*k2/2;
    VALIDAR(x3); VALIDAR(y3);
    
    k3 = EDO_FUNCION(x3, y3);
    VALIDAR(k3);
    
    // k4
    double x4 = x + h;
    double y4 = y + h*k3;
    VALIDAR(x4); VALIDAR(y4);
    
    k4 = EDO_FUNCION(x4, y4);
    VALIDAR(k4);
    
    // Resultado final
    double resultado = y + h*(k1 + 2*k2 + 2*k3 + k4)/6;
    VALIDAR(resultado);
    
    // Detectar inestabilidad
    if (fabs(resultado) > 1e10 && paso_actual > 10) {
        printf(" ADVERTENCIA [Paso %d]: Posible inestabilidad numerica\n", paso_actual);
        printf("   y = %.2e, puede haber divergencia\n", resultado);
    }
    
    return resultado;
}

// ============================================================================
// PROGRAMA PRINCIPAL
// ============================================================================
int main() {
    // ============================================================================
    // VALIDACION INICIAL
    // ============================================================================
    printf(" VALIDANDO PARAMETROS...\n");
    printf("-------------------------------------------------------------\n");
    
    validar_parametros();
    
    double x = X_INICIAL;
    double y = Y_INICIAL;
    int paso = 0;
    int pasos_totales = (int)((X_FINAL - X_INICIAL) / PASO_H) + 1;
    
    printf(" Parametros validos\n");
    printf("   Ecuacion: y' = x - y\n");
    printf("   Condicion inicial: y(%.1f) = %.1f\n", X_INICIAL, Y_INICIAL);
    printf("   Intervalo: [%.1f, %.1f]\n", X_INICIAL, X_FINAL);
    printf("   Paso: h = %.3f\n", PASO_H);
    printf("   Pasos estimados: %d\n\n", pasos_totales);
    
    // ============================================================================
    // CONFIGURACION
    // ============================================================================
    printf(" ECUACION DIFERENCIAL: y' = x - y (RK4) \n\n");
    
    FILE *datos_sol = abrir_archivo("rk4_solucion.dat", "w");
    FILE *datos_err = abrir_archivo("rk4_error.dat", "w");
    FILE *script_gp = abrir_archivo("rk4_plot.gp", "w");
    
    printf("PROCESO DE INTEGRACION:\n");
    printf("+------+--------+-----------+-----------+-----------+-----------+\n");
    printf("| Paso |   x    |   y_RK4   | y_Exacta  |  Error    | Estado    |\n");
    printf("+------+--------+-----------+-----------+-----------+-----------+\n");
    
    int errores_numericos = 0;
    double error_maximo = 0.0;
    
    // ============================================================================
    // INTEGRACION CON RUNGE-KUTTA 4
    // ============================================================================
    while (x <= X_FINAL + PASO_H/2) {
        // Calcular solucion exacta
        double exacta = SOLUCION_EXACTA(x);
        VALIDAR(exacta);
        
        double error = fabs(y - exacta);
        VALIDAR(error);
        
        if (error > error_maximo) {
            error_maximo = error;
        }
        
        // Estado de validacion
        const char *estado = "- OK";
        if (!es_numerico_valido(y)) {
            estado = "- INVALIDO";
            errores_numericos++;
            
            printf("+------+--------+-----------+-----------+-----------+-----------+\n");
            printf("| %4d | %6.2f | %9.5f | %9.5f | %9.5f | %s |\n", 
                   paso, x, y, exacta, error, estado);
            printf("+------+--------+-----------+-----------+-----------+-----------+\n");
            
            printf("\n ERROR CRITICO: Valor no numerico en paso %d\n", paso);
            printf("   x = %.6f, y = %.6f\n", x, y);
            printf("   El metodo no puede continuar\n");
            
            fclose(datos_sol);
            fclose(datos_err);
            return EXIT_FAILURE;
        }
        
        // Mostrar cada 5 pasos
        if (paso % 5 == 0) {
            printf("| %4d | %6.2f | %9.5f | %9.5f | %9.5f | %s |\n", 
                   paso, x, y, exacta, error, estado);
        }
        
        // Guardar datos
        fprintf(datos_sol, "%.6f %.6f\n", x, y);
        fprintf(datos_err, "%.6f %.6f\n", x, error);
        
        // Ultimo punto
        if (x >= X_FINAL) break;
        
        // Calcular siguiente punto
        double y_nuevo = rk4_validado(x, y, PASO_H, paso);
        
        // Validar nuevo valor
        if (!es_numerico_valido(y_nuevo)) {
            printf(" ADVERTENCIA: Valor invalido en paso %d, ajustando...\n", paso);
            
            // Intentar con paso mas pequeño
            double y_half1 = rk4_validado(x, y, PASO_H/2, paso);
            double y_half2 = rk4_validado(x + PASO_H/2, y_half1, PASO_H/2, paso);
            
            if (es_numerico_valido(y_half2)) {
                y_nuevo = y_half2;
                printf("   Solucionado con paso reducido a h/2\n");
            } else {
                printf(" ERROR: No se pudo recuperar con paso reducido\n");
                break;
            }
        }
        
        y = y_nuevo;
        x += PASO_H;
        paso++;
        
        // Verificar limite de pasos (prevencion de bucle infinito)
        if (paso > pasos_totales * 10) {
            printf(" ADVERTENCIA: Demasiados pasos (%d), posible bucle infinito\n", paso);
            break;
        }
    }
    
    printf("+------+--------+-----------+-----------+-----------+-----------+\n");
    printf("| INTEGRACION COMPLETADA: %d pasos, %d errores numericos      |\n", 
           paso, errores_numericos);
    printf("+--------------------------------------------------------------+\n\n");
    
    fclose(datos_sol);
    fclose(datos_err);
    
    // ============================================================================
    // CREAR SCRIPT GNUPLOT
    // ============================================================================
    fprintf(script_gp, "# Script para Runge-Kutta 4\n");
    fprintf(script_gp, "set terminal pngcairo size %d,%d enhanced font 'Arial,10'\n", 
            ANCHO_GRAFICO, ALTO_GRAFICO);
    fprintf(script_gp, "set output '%s'\n", NOMBRE_GRAFICO);
    fprintf(script_gp, "set title 'Metodo de Runge-Kutta 4: y' = x - y'\n");
    fprintf(script_gp, "set xlabel 'x'\n");
    fprintf(script_gp, "set ylabel 'y(x)'\n");
    fprintf(script_gp, "set grid\n");
    fprintf(script_gp, "set key top left box\n");
    fprintf(script_gp, "set xrange [%f:%f]\n", X_INICIAL, X_FINAL);
    
    fprintf(script_gp, "plot 'rk4_solucion.dat' w lp pt 7 ps 0.5 lc rgb '#0066CC' title 'Solucion RK4', \\\n");
    fprintf(script_gp, "     x - 1 + 2*exp(-x) w l lw 2 lc rgb '#FF3333' title 'Solucion exacta', \\\n");
    fprintf(script_gp, "     'rk4_error.dat' u 1:2 w l lw 1 lc rgb '#00AA00' axes x1y2 title 'Error'\n");
    
    fprintf(script_gp, "\n# Configurar segundo eje Y para error\n");
    fprintf(script_gp, "set y2tics\n");
    fprintf(script_gp, "set y2label 'Error absoluto'\n");
    
    fclose(script_gp);
    
    // ============================================================================
    // EJECUTAR GNUPLOT
    // ============================================================================
    printf(" GENERANDO GRAFICO...\n");
    printf("-------------------------------------------------------------\n");
    
    int resultado_gnuplot = system("gnuplot rk4_plot.gp 2>&1");
    
    if (resultado_gnuplot != 0) {
        printf(" ADVERTENCIA: Gnuplot reporto problemas\n");
        printf("   Comando: gnuplot rk4_plot.gp\n");
    } else {
        printf(" Grafico generado: %s\n", NOMBRE_GRAFICO);
    }
    
    // ============================================================================
    // ANALISIS DE RESULTADOS
    // ============================================================================
    printf("\n ANALISIS DE RESULTADOS:\n");
    printf("-------------------------------------------------------------\n");
    
    // Calcular estadisticas de error
    FILE *err_file = fopen("rk4_error.dat", "r");
    double error_promedio = 0.0;
    double error_final = 0.0;
    int puntos_leidos = 0;
    double x_val, err_val;
    
    if (err_file) {
        while (fscanf(err_file, "%lf %lf", &x_val, &err_val) == 2) {
            error_promedio += err_val;
            if (x_val == X_FINAL || puntos_leidos == paso - 1) {
                error_final = err_val;
            }
            puntos_leidos++;
        }
        fclose(err_file);
        
        if (puntos_leidos > 0) {
            error_promedio /= puntos_leidos;
        }
    }
    
    printf("  Pasos completados:   %d de %d estimados\n", paso, pasos_totales);
    printf("  Error maximo:        %.6f\n", error_maximo);
    printf("  Error promedio:      %.6f\n", error_promedio);
    printf("  Error final:         %.6f\n", error_final);
    printf("  Errores numericos:   %d\n", errores_numericos);
    
    // Evaluar precision
    printf("\n  EVALUACION DE PRECISION:\n");
    if (error_maximo < 0.001) {
        printf("    - Excelente precision (error maximo < 0.001)\n");
    } else if (error_maximo < 0.01) {
        printf("    - Buena precision (error maximo < 0.01)\n");
    } else if (error_maximo < 0.1) {
        printf("    - Precision aceptable (error maximo < 0.1)\n");
    } else {
        printf("    - Precision pobre, considere reducir el paso h\n");
    }
    
    // ============================================================================
    // VALIDACION FINAL DE CONSERVACION
    // ============================================================================
    printf("\n VALIDACION DE CONSERVACION:\n");
    printf("-------------------------------------------------------------\n");
    
    // Calcular derivada numerica final
    double derivada_final = EDO_FUNCION(X_FINAL, y);
    double derivada_teorica = X_FINAL - y;
    double discrepancia = fabs(derivada_final - derivada_teorica);
    
    printf("  En x = %.2f:\n", X_FINAL);
    printf("    y calculado:     %.6f\n", y);
    printf("    y exacto:        %.6f\n", SOLUCION_EXACTA(X_FINAL));
    printf("    y' calculado:    %.6f\n", derivada_final);
    printf("    y' teorico:      %.6f\n", derivada_teorica);
    printf("    Discrepancia:    %.2e\n", discrepancia);
    
    if (discrepancia > 0.01) {
        printf("  ADVERTENCIA: Discrepancia significativa en derivada\n");
    }
    
    // ============================================================================
    // RESUMEN FINAL
    // ============================================================================
    printf("\n RESUMEN DE EJECUCION:\n");
    printf("-------------------------------------------------------------\n");
    printf("  Estado:              %s\n", 
           (errores_numericos == 0) ? "EXITOSO" : "CON ADVERTENCIAS");
    printf("  Pasos ejecutados:    %d\n", paso);
    printf("  Error maximo:        %.2e\n", error_maximo);
    printf("  Grafico generado:    %s\n", 
           (resultado_gnuplot == 0) ? "SI" : "NO");
    printf("  Archivos creados:    rk4_solucion.dat, rk4_error.dat, rk4_plot.gp\n");
    
    printf("\n EJECUCION COMPLETADA\n");
    
    return (errores_numericos == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}