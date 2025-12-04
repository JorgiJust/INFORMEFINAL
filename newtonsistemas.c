// 2_newton_sistemas.c
// Newton para sistemas con validaciones robustas

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <errno.h>

// ============================================================================
// PARAMETROS CONFIGURABLES
// ============================================================================
#define F1(x,y)             ((x)*(x) + (y)*(y) - 4)
#define F2(x,y)             (exp(x) + (y) - 1)
#define DF1_DX(x,y)         (2*(x))
#define DF1_DY(x,y)         (2*(y))
#define DF2_DX(x,y)         (exp(x))
#define DF2_DY(x,y)         (1)
#define X_INICIAL           1.0
#define Y_INICIAL           1.0
#define TOLERANCIA          1e-6
#define MAX_ITER            50
#define GRAFICO_RANGO_X     3.0
#define GRAFICO_RANGO_Y     3.0
#define GRAFICO_PUNTOS      200
#define NOMBRE_GRAFICO      "sistema_grafico.png"
#define ANCHO_GRAFICO       900
#define ALTO_GRAFICO        700
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
    if (!es_numerico_valido(valor)) {
        printf("ERROR en linea %d: %s = Valor invalido\n", linea, nombre);
        exit(EXIT_FAILURE);
    }
}

#define VALIDAR(variable) verificar_nan_inf(#variable, variable, __LINE__)

void validar_punto(double x, double y, const char *contexto) {
    if (!es_numerico_valido(x) || !es_numerico_valido(y)) {
        printf("ERROR en %s: Punto invalido (%.6f, %.6f)\n", contexto, x, y);
        exit(EXIT_FAILURE);
    }
}

FILE* abrir_archivo(const char *nombre, const char *modo) {
    FILE *archivo = fopen(nombre, modo);
    if (archivo == NULL) {
        printf("ERROR: No se pudo abrir '%s' (modo: %s)\n", nombre, modo);
        printf("   errno: %d\n", errno);
        exit(EXIT_FAILURE);
    }
    return archivo;
}

// ============================================================================
// FUNCIONES PRINCIPALES
// ============================================================================
void generar_datos_curvas() {
    FILE *curvas = abrir_archivo("sistema_curvas.dat", "w");
    
    fprintf(curvas, "# Curva 1: x^2 + y^2 = 4\n");
    for (int i = 0; i <= GRAFICO_PUNTOS; i++) {
        double t = 2 * M_PI * i / GRAFICO_PUNTOS;
        double x = 2 * cos(t);
        double y = 2 * sin(t);
        VALIDAR(x); VALIDAR(y);
        fprintf(curvas, "%.6f %.6f\n", x, y);
    }
    fprintf(curvas, "\n\n# Curva 2: e^x + y = 1\n");
    for (int i = 0; i <= GRAFICO_PUNTOS; i++) {
        double xi = -GRAFICO_RANGO_X + (2*GRAFICO_RANGO_X * i / GRAFICO_PUNTOS);
        double yi = 1 - exp(xi);
        VALIDAR(xi); VALIDAR(yi);
        fprintf(curvas, "%.6f %.6f\n", xi, yi);
    }
    fclose(curvas);
}

void crear_script_gnuplot(double sol_x, double sol_y) {
    FILE *gp = abrir_archivo("sistema_plot.gp", "w");
    
    fprintf(gp, "# Script para sistema de ecuaciones\n");
    fprintf(gp, "set terminal pngcairo size %d,%d enhanced font 'Arial,10'\n", 
            ANCHO_GRAFICO, ALTO_GRAFICO);
    fprintf(gp, "set output '%s'\n", NOMBRE_GRAFICO);
    fprintf(gp, "set title 'Sistema: x^2+y^2=4 y e^x+y=1'\n");
    fprintf(gp, "set xlabel 'x'\n");
    fprintf(gp, "set ylabel 'y'\n");
    fprintf(gp, "set grid\n");
    fprintf(gp, "set size ratio -1\n");
    fprintf(gp, "set xrange [%lf:%lf]\n", -GRAFICO_RANGO_X, GRAFICO_RANGO_X);
    fprintf(gp, "set yrange [%lf:%lf]\n", -GRAFICO_RANGO_Y, GRAFICO_RANGO_Y);
    fprintf(gp, "set key box opaque\n\n");
    
    fprintf(gp, "plot 'sistema_curvas.dat' index 0 w l lw 2 lc rgb '#0066CC' title 'x^2 + y^2 = 4', \\\n");
    fprintf(gp, "     'sistema_curvas.dat' index 1 w l lw 2 lc rgb '#CC0066' title 'e^x + y = 1', \\\n");
    fprintf(gp, "     'sistema_trayectoria.dat' w l lw 1.5 lc rgb '#00AA00' title 'Trayectoria Newton', \\\n");
    fprintf(gp, "     'sistema_trayectoria.dat' w p pt 7 ps 1 lc rgb '#00AA00' notitle, \\\n");
    fprintf(gp, "     %lf, %lf w p pt 9 ps 2 lc rgb '#000000' title 'Solucion: (%.4f, %.4f)'\n", 
            sol_x, sol_y, sol_x, sol_y);
    
    fclose(gp);
}

int ejecutar_gnuplot() {
    printf("\nGenerando grafico...\n");
    int resultado = system("gnuplot sistema_plot.gp 2>&1");
    
    if (resultado != 0) {
        printf("ADVERTENCIA: Error al ejecutar Gnuplot\n");
        printf("   Comando fallido: gnuplot sistema_plot.gp\n");
        return 0;
    }
    
    printf("EXITO: Grafico generado: %s\n", NOMBRE_GRAFICO);
    return 1;
}

int main() {
    double x = X_INICIAL, y = Y_INICIAL, error;
    int iteracion = 0;
    
    // ============================================================================
    // VALIDACION INICIAL
    // ============================================================================
    printf("Validando parametros iniciales...\n");
    
    validar_punto(x, y, "punto inicial");
    
    double f1_inicial = F1(x, y);
    double f2_inicial = F2(x, y);
    VALIDAR(f1_inicial);
    VALIDAR(f2_inicial);
    
    if (TOLERANCIA <= 0) {
        printf("ERROR: TOLERANCIA debe ser positiva\n");
        return EXIT_FAILURE;
    }
    
    printf("EXITO: Validacion inicial exitosa\n");
    printf("   f1(%.1f, %.1f) = %.3f\n", x, y, f1_inicial);
    printf("   f2(%.1f, %.1f) = %.3f\n\n", x, y, f2_inicial);
    
    // ============================================================================
    // CONFIGURACION
    // ============================================================================
    printf("===============================================================\n");
    printf("          SISTEMA DE ECUACIONES NO LINEALES (2D)              \n");
    printf("===============================================================\n\n");
    
    printf("PROCESO DE CALCULO:\n");
    printf("================================================================================\n");
    printf("| Iter|     x     |     y     |    f1     |    f2     |  Det(J)   |  Error    |\n");
    printf("================================================================================\n");
    
    // Archivos
    FILE *datos_iter = abrir_archivo("sistema_iteraciones.dat", "w");
    FILE *datos_tray = abrir_archivo("sistema_trayectoria.dat", "w");
    
    fprintf(datos_iter, "# iter x y f1 f2 det_j error\n");
    fprintf(datos_tray, "%.6f %.6f\n", x, y);
    
    // ============================================================================
    // METODO DE NEWTON CON VALIDACIONES
    // ============================================================================
    do {
        double f1 = F1(x, y);
        double f2 = F2(x, y);
        VALIDAR(f1); VALIDAR(f2);
        
        // Jacobiano
        double df1_dx = DF1_DX(x, y);
        double df1_dy = DF1_DY(x, y);
        double df2_dx = DF2_DX(x, y);
        double df2_dy = DF2_DY(x, y);
        
        VALIDAR(df1_dx); VALIDAR(df1_dy);
        VALIDAR(df2_dx); VALIDAR(df2_dy);
        
        double det = df1_dx*df2_dy - df1_dy*df2_dx;
        VALIDAR(det);
        
        // Validacion de Jacobiano
        if (fabs(det) < 1e-15) {
            printf("================================================================================\n");
            printf("| ERROR CRITICO: Jacobiano singular                                          |\n");
            printf("|   det(J) = %.2e en (%.6f, %.6f)                             |\n", det, x, y);
            printf("|   f1 = %.6f, f2 = %.6f                                       |\n", f1, f2);
            printf("================================================================================\n");
            fclose(datos_iter);
            fclose(datos_tray);
            return EXIT_FAILURE;
        }
        
        // Resolver sistema
        double dx = (-f1*df2_dy + f2*df1_dy) / det;
        double dy = (-df1_dx*f2 + f1*df2_dx) / det;
        
        VALIDAR(dx); VALIDAR(dy);
        
        error = sqrt(dx*dx + dy*dy);
        VALIDAR(error);
        
        // Mostrar resultados
        printf("| %3d | %9.6f | %9.6f | %9.6f | %9.6f | %9.2e | %9.6f |\n", 
               iteracion, x, y, f1, f2, det, error);
        
        // Guardar
        fprintf(datos_iter, "%d %.6f %.6f %.6f %.6f %.6e %.6f\n", 
                iteracion, x, y, f1, f2, det, error);
        
        // Actualizar con validacion
        double x_nuevo = x + dx;
        double y_nuevo = y + dy;
        
        validar_punto(x_nuevo, y_nuevo, "nuevo punto");
        
        x = x_nuevo;
        y = y_nuevo;
        
        fprintf(datos_tray, "%.6f %.6f\n", x, y);
        iteracion++;
        
        // Deteccion de divergencia
        if (error > 1e5 && iteracion > 3) {
            printf("================================================================================\n");
            printf("| ADVERTENCIA: Posible divergencia                                           |\n");
            printf("|   Error creciente: %.2e                                                    |\n", error);
            printf("================================================================================\n");
            break;
        }
        
        // Verificar convergencia
        if (error < TOLERANCIA) {
            printf("================================================================================\n");
            printf("| EXITO: CONVERGENCIA ALCANZADA                                              |\n");
            printf("|   Error: %.2e < Tolerancia: %.2e                                       |\n", 
                   error, TOLERANCIA);
            printf("================================================================================\n");
            break;
        }
        
        if (iteracion >= MAX_ITER) {
            printf("================================================================================\n");
            printf("| ADVERTENCIA: LIMITE DE ITERACIONES                                         |\n");
            printf("|   No convergio en %d iteraciones                                          |\n", MAX_ITER);
            printf("|   Ultimo error: %.2e                                                    |\n", error);
            printf("================================================================================\n");
            break;
        }
        
    } while (1);
    
    fclose(datos_iter);
    fclose(datos_tray);
    
    // ============================================================================
    // VALIDACION DE SOLUCION FINAL
    // ============================================================================
    printf("\nValidando solucion final...\n");
    
    double f1_final = F1(x, y);
    double f2_final = F2(x, y);
    VALIDAR(f1_final); VALIDAR(f2_final);
    
    double error_f1 = fabs(f1_final);
    double error_f2 = fabs(f2_final);
    
    if (error_f1 > 0.01 || error_f2 > 0.01) {
        printf("ADVERTENCIA: La solucion no satisface bien las ecuaciones\n");
        printf("   f1(x,y) = %.2e (deberia ser ~0)\n", error_f1);
        printf("   f2(x,y) = %.2e (deberia ser ~0)\n", error_f2);
    } else {
        printf("EXITO: Solucion valida las ecuaciones\n");
    }
    
    // ============================================================================
    // GENERAR GRAFICOS
    // ============================================================================
    generar_datos_curvas();
    crear_script_gnuplot(x, y);
    int grafico_ok = ejecutar_gnuplot();
    
    // ============================================================================
    // RESULTADOS FINALES
    // ============================================================================
    printf("\nRESULTADOS FINALES:\n");
    printf("-----------------------------------------------------------------\n");
    printf("  Solucion:         x = %.8f, y = %.8f\n", x, y);
    printf("  f1(x,y) =         %.2e\n", f1_final);
    printf("  f2(x,y) =         %.2e\n", f2_final);
    printf("  Iteraciones:      %d de %d\n", iteracion, MAX_ITER);
    printf("  Error final:      %.2e\n", error);
    printf("  Estado:           %s\n", 
           (error < TOLERANCIA) ? "CONVERGENCIA" : "ITERACIONES MAXIMAS");
    
    printf("\nARCHIVOS GENERADOS:\n");
    printf("-----------------------------------------------------------------\n");
    printf("  EXITO: sistema_iteraciones.dat -> %d iteraciones\n", iteracion);
    printf("  EXITO: sistema_trayectoria.dat -> Trayectoria completa\n");
    printf("  EXITO: sistema_curvas.dat      -> Curvas de ecuaciones\n");
    printf("  EXITO: sistema_plot.gp         -> Script Gnuplot\n");
    if (grafico_ok) {
        printf("  EXITO: %s       -> Grafico final\n", NOMBRE_GRAFICO);
    }
    
    printf("\n===============================================================\n");
    printf("                      EJECUCION COMPLETADA                     \n");
    printf("===============================================================\n");
    
    return EXIT_SUCCESS;
}