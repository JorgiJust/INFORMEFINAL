// 6_ecuacion2.c
// EDO y'' + y = 0 con Runge-Kutta 4 y validaciones robustas

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <errno.h>

// ============================================================================
// PARAMETROS CONFIGURABLES
// ============================================================================
#define EDO_FUNCION(x,y,yp) (-(y))
#define SOLUCION_EXACTA(x)  (sin(x))
#define X_INICIAL           0.0
#define X_FINAL             4*M_PI
#define Y_INICIAL           0.0
#define YP_INICIAL          1.0
#define PASO_H              0.05
#define NOMBRE_GRAFICO      "ypp_grafico.png"
#define ANCHO_GRAFICO       800
#define ALTO_GRAFICO        1000
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
        printf("   Posible causa: Division por cero o operacion invalida\n");
        exit(EXIT_FAILURE);
    }
    if (isinf(valor)) {
        printf(" ERROR [Linea %d]: %s = Infinito\n", linea, nombre);
        printf("   Posible causa: Overflow numerico\n");
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
        exit(EXIT_FAILURE);
    }
    
    if (!es_numerico_valido(Y_INICIAL) || !es_numerico_valido(YP_INICIAL)) {
        printf(" ERROR: Condiciones iniciales invalidas\n");
        printf("   y(0) = %f, y'(0) = %f\n", Y_INICIAL, YP_INICIAL);
        exit(EXIT_FAILURE);
    }
    
    // Validar periodicidad esperada
    if (fabs(X_FINAL - 4*M_PI) > 1e-10) {
        printf(" ADVERTENCIA: X_FINAL no es multiplo exacto de π\n");
        printf("   Para solucion periodica, use X_FINAL = n*π\n");
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
// RUNGE-KUTTA 4 PARA SISTEMAS CON VALIDACION
// ============================================================================
void rk4_sistema_validado(double x, double *y, double *yp, double h, int paso_actual) {
    double k1_y, k1_yp, k2_y, k2_yp, k3_y, k3_yp, k4_y, k4_yp;
    
    // k1
    k1_y = *yp;
    k1_yp = EDO_FUNCION(x, *y, *yp);
    VALIDAR(k1_y); VALIDAR(k1_yp);
    
    // k2
    double x2 = x + h/2;
    double y2 = *y + h*k1_y/2;
    double yp2 = *yp + h*k1_yp/2;
    VALIDAR(x2); VALIDAR(y2); VALIDAR(yp2);
    
    k2_y = yp2;
    k2_yp = EDO_FUNCION(x2, y2, yp2);
    VALIDAR(k2_y); VALIDAR(k2_yp);
    
    // k3
    double x3 = x + h/2;
    double y3 = *y + h*k2_y/2;
    double yp3 = *yp + h*k2_yp/2;
    VALIDAR(x3); VALIDAR(y3); VALIDAR(yp3);
    
    k3_y = yp3;
    k3_yp = EDO_FUNCION(x3, y3, yp3);
    VALIDAR(k3_y); VALIDAR(k3_yp);
    
    // k4
    double x4 = x + h;
    double y4 = *y + h*k3_y;
    double yp4 = *yp + h*k3_yp;
    VALIDAR(x4); VALIDAR(y4); VALIDAR(yp4);
    
    k4_y = yp4;
    k4_yp = EDO_FUNCION(x4, y4, yp4);
    VALIDAR(k4_y); VALIDAR(k4_yp);
    
    // Actualizar valores
    double y_nuevo = *y + h*(k1_y + 2*k2_y + 2*k3_y + k4_y)/6;
    double yp_nuevo = *yp + h*(k1_yp + 2*k2_yp + 2*k3_yp + k4_yp)/6;
    
    VALIDAR(y_nuevo); VALIDAR(yp_nuevo);
    
    // Verificar conservacion de energia (E = y² + y'²)
    double energia_antes = (*y)*(*y) + (*yp)*(*yp);
    double energia_despues = y_nuevo*y_nuevo + yp_nuevo*yp_nuevo;
    double delta_energia = fabs(energia_despues - energia_antes);
    
    if (delta_energia > 0.01 && paso_actual > 10) {
        printf(" ADVERTENCIA [Paso %d]: Energia no se conserva bien\n", paso_actual);
        printf("   ΔE = %.2e (deberia ser ~0)\n", delta_energia);
    }
    
    // Verificar crecimiento exponencial (indicador de inestabilidad)
    if (fabs(y_nuevo) > 10*fabs(*y) && paso_actual > 5) {
        printf(" ADVERTENCIA [Paso %d]: Posible inestabilidad\n", paso_actual);
        printf("   y crecio de %.2e a %.2e\n", *y, y_nuevo);
    }
    
    *y = y_nuevo;
    *yp = yp_nuevo;
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
    double yp = YP_INICIAL;
    int paso = 0;
    int pasos_totales = (int)((X_FINAL - X_INICIAL) / PASO_H) + 1;
    
    printf(" Parametros validos\n");
    printf("   Ecuacion: y'' + y = 0\n");
    printf("   Condiciones: y(0) = %.1f, y'(0) = %.1f\n", Y_INICIAL, YP_INICIAL);
    printf("   Intervalo: [%.1f, %.1f]\n", X_INICIAL, X_FINAL);
    printf("   Paso: h = %.3f\n", PASO_H);
    printf("   Pasos estimados: %d\n\n", pasos_totales);
    
    // ============================================================================
    // CONFIGURACION
    // ============================================================================
    printf(" ECUACION DIFERENCIAL: y'' + y = 0 (RK4) \n\n");
    
    FILE *datos_sol = abrir_archivo("ypp_solucion.dat", "w");
    FILE *datos_der = abrir_archivo("ypp_derivada.dat", "w");
    FILE *datos_fase = abrir_archivo("ypp_fase.dat", "w");
    FILE *script_gp = abrir_archivo("ypp_plot.gp", "w");
    
    printf("PROCESO DE INTEGRACION:\n");
    printf("+------+--------+-----------+-----------+-----------+-----------+\n");
    printf("| Paso |   x    |   y(x)    |   y'(x)   |  Error    |  Energia  |\n");
    printf("+------+--------+-----------+-----------+-----------+-----------+\n");
    
    int errores_numericos = 0;
    double error_maximo = 0.0;
    double energia_inicial = Y_INICIAL*Y_INICIAL + YP_INICIAL*YP_INICIAL;
    
    // ============================================================================
    // INTEGRACION CON RUNGE-KUTTA 4
    // ============================================================================
    while (x <= X_FINAL + PASO_H/2) {
        // Calcular solucion exacta y error
        double exacta = SOLUCION_EXACTA(x);
        VALIDAR(exacta);
        
        double error = fabs(y - exacta);
        VALIDAR(error);
        
        if (error > error_maximo) {
            error_maximo = error;
        }
        
        // Calcular energia actual
        double energia_actual = y*y + yp*yp;
        VALIDAR(energia_actual);
        
        // Estado de validacion
        const char *estado = "- OK";
        if (!es_numerico_valido(y) || !es_numerico_valido(yp)) {
            estado = "- INVALIDO";
            errores_numericos++;
            
            printf("+------+--------+-----------+-----------+-----------+-----------+\n");
            printf("| %4d | %6.3f | %9.5f | %9.5f | %9.5f | %9.5f | %s |\n", 
                   paso, x, y, yp, error, energia_actual, estado);
            printf("+------+--------+-----------+-----------+-----------+-----------+\n");
            
            printf("\n ERROR CRITICO: Valores no numericos en paso %d\n", paso);
            printf("   x = %.6f, y = %.6f, y' = %.6f\n", x, y, yp);
            
            fclose(datos_sol);
            fclose(datos_der);
            fclose(datos_fase);
            return EXIT_FAILURE;
        }
        
        // Mostrar cada 20 pasos
        if (paso % 20 == 0) {
            printf("| %4d | %6.3f | %9.5f | %9.5f | %9.5f | %9.5f | %s |\n", 
                   paso, x, y, yp, error, energia_actual, estado);
        }
        
        // Guardar datos
        fprintf(datos_sol, "%.6f %.6f\n", x, y);
        fprintf(datos_der, "%.6f %.6f\n", x, yp);
        fprintf(datos_fase, "%.6f %.6f\n", y, yp);
        
        // Ultimo punto
        if (x >= X_FINAL) break;
        
        // Calcular siguiente punto
        rk4_sistema_validado(x, &y, &yp, PASO_H, paso);
        
        x += PASO_H;
        paso++;
        
        // Verificar limite de pasos
        if (paso > pasos_totales * 10) {
            printf(" ADVERTENCIA: Demasiados pasos (%d)\n", paso);
            break;
        }
    }
    
    printf("+------+--------+-----------+-----------+-----------+-----------+\n");
    printf("| INTEGRACION COMPLETADA: %d pasos                            |\n", paso);
    printf("+--------------------------------------------------------------+\n\n");
    
    fclose(datos_sol);
    fclose(datos_der);
    fclose(datos_fase);
    
    // ============================================================================
    // CREAR SCRIPT GNUPLOT
    // ============================================================================
    fprintf(script_gp, "# Script para ecuacion y'' + y = 0\n");
    fprintf(script_gp, "set terminal pngcairo size %d,%d enhanced font 'Arial,10'\n", 
            ANCHO_GRAFICO, ALTO_GRAFICO);
    fprintf(script_gp, "set output '%s'\n", NOMBRE_GRAFICO);
    
    fprintf(script_gp, "\n# Configurar multiples graficos\n");
    fprintf(script_gp, "set multiplot layout 2,1\n");
    fprintf(script_gp, "set lmargin 10\n");
    fprintf(script_gp, "set rmargin 5\n\n");
    
    // Grafico 1: Solucion
    fprintf(script_gp, "# Grafico 1: Solucion y(x)\n");
    fprintf(script_gp, "set title 'Solucion: y'' + y = 0'\n");
    fprintf(script_gp, "set xlabel 'x'\n");
    fprintf(script_gp, "set ylabel 'y(x)'\n");
    fprintf(script_gp, "set grid\n");
    fprintf(script_gp, "set key top left box\n");
    fprintf(script_gp, "plot 'ypp_solucion.dat' w l lw 2 lc rgb '#0066CC' title 'Solucion RK4', \\\n");
    fprintf(script_gp, "     sin(x) w l lw 2 lc rgb '#FF3333' dt 2 title 'sin(x) (exacta)'\n\n");
    
    // Grafico 2: Plano de fase
    fprintf(script_gp, "# Grafico 2: Plano de fase\n");
    fprintf(script_gp, "set title 'Plano de fase: y vs y''\n");
    fprintf(script_gp, "set xlabel 'y(x)'\n");
    fprintf(script_gp, "set ylabel 'y'(x)'\n");
    fprintf(script_gp, "set grid\n");
    fprintf(script_gp, "set key off\n");
    fprintf(script_gp, "set size ratio -1\n");
    fprintf(script_gp, "plot 'ypp_fase.dat' w l lw 1.5 lc rgb '#00AA00' title 'Trayectoria'\n\n");
    
    fprintf(script_gp, "unset multiplot\n");
    fclose(script_gp);
    
    // ============================================================================
    // EJECUTAR GNUPLOT
    // ============================================================================
    printf(" GENERANDO GRAFICO...\n");
    printf("-------------------------------------------------------------\n");
    
    int resultado_gnuplot = system("gnuplot ypp_plot.gp 2>&1");
    
    if (resultado_gnuplot != 0) {
        printf(" ADVERTENCIA: Gnuplot reporto problemas\n");
    } else {
        printf(" Grafico generado: %s\n", NOMBRE_GRAFICO);
    }
    
    // ============================================================================
    // ANALISIS DE RESULTADOS
    // ============================================================================
    printf("\n ANALISIS DE RESULTADOS:\n");
    printf("-------------------------------------------------------------\n");
    
    // Calcular energia final y variacion
    double energia_final = y*y + yp*yp;
    double variacion_energia = fabs(energia_final - energia_inicial);
    double variacion_relativa = 100 * variacion_energia / energia_inicial;
    
    // Calcular periodicidad
    double periodo_teorico = 2*M_PI;
    int ciclos_completos = (int)(X_FINAL / periodo_teorico);
    
    printf("  Pasos completados:   %d\n", paso);
    printf("  Error maximo:        %.6f\n", error_maximo);
    printf("  Energia inicial:     %.6f\n", energia_inicial);
    printf("  Energia final:       %.6f\n", energia_final);
    printf("  Variacion energia:   %.2e (%.2f%%)\n", 
           variacion_energia, variacion_relativa);
    printf("  Ciclos completos:    %d\n", ciclos_completos);
    printf("  Errores numericos:   %d\n", errores_numericos);
    
    // Evaluar conservacion de energia
    printf("\n  EVALUACION DE CONSERVACION DE ENERGIA:\n");
    if (variacion_relativa < 0.1) {
        printf("    - Excelente conservacion (ΔE < 0.1%%)\n");
    } else if (variacion_relativa < 1.0) {
        printf("    - Buena conservacion (ΔE < 1%%)\n");
    } else if (variacion_relativa < 5.0) {
        printf("    - Conservacion aceptable (ΔE < 5%%)\n");
    } else {
        printf("    - Mala conservacion, metodo puede ser inestable\n");
    }
    
    // Evaluar periodicidad
    printf("\n  EVALUACION DE PERIODICIDAD:\n");
    double y_final_teorico = sin(X_FINAL);
    double error_periodicidad = fabs(y - y_final_teorico);
    
    if (error_periodicidad < 0.01) {
        printf("    - Buena periodicidad (error < 0.01)\n");
    } else if (error_periodicidad < 0.1) {
        printf("    - Periodicidad aceptable (error < 0.1)\n");
    } else {
        printf("    - Periodicidad pobre, posible acumulacion de error\n");
    }
    
    // ============================================================================
    // VALIDACION DE PROPIEDADES MATEMATICAS
    // ============================================================================
    printf("\n VALIDACION DE PROPIEDADES MATEMATICAS:\n");
    printf("-------------------------------------------------------------\n");
    
    // Verificar que satisface la EDO
    double ypp_numerica = EDO_FUNCION(x, y, yp);
    double residual = ypp_numerica + y;  // y'' + y deberia ser 0
    
    printf("  En x = %.4f:\n", x);
    printf("    y calculado:       %.8f\n", y);
    printf("    y' calculado:      %.8f\n", yp);
    printf("    y'' calculado:     %.8f\n", ypp_numerica);
    printf("    Residual (y''+y):  %.2e (deberia ser ~0)\n", residual);
    
    if (fabs(residual) > 0.1) {
        printf("  ADVERTENCIA: Residual grande, solucion puede no satisfacer EDO\n");
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
    printf("  ΔEnergia:            %.2e (%.2f%%)\n", 
           variacion_energia, variacion_relativa);
    printf("  Grafico generado:    %s\n", 
           (resultado_gnuplot == 0) ? "SI" : "NO");
    printf("  Archivos creados:    ypp_solucion.dat, ypp_derivada.dat, ypp_fase.dat\n");
    
    printf("\n EJECUCION COMPLETADA\n");
    
    return (errores_numericos == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}