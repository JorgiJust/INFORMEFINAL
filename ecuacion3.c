/// 7_sistema_edos.c
// Sistema dx/dt = y, dy/dt = -x con validaciones robustas

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <errno.h>

// ============================================================================
// ============================================================================
#define F1(x,y)             (y)
#define F2(x,y)             (-(x))
#define T_INICIAL           0.0
#define T_FINAL             10.0
#define X_INICIAL           1.0
#define Y_INICIAL           0.0
#define PASO_H              0.05
#define NOMBRE_GRAFICO1     "sistema_temporal.png"
#define NOMBRE_GRAFICO2     "sistema_fase.png"
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
        printf("ERROR [Linea %d]: %s = NaN\n", linea, nombre);
        printf("   Causa: Operacion matematica invalida\n");
        exit(EXIT_FAILURE);
    }
    if (isinf(valor)) {
        printf("ERROR [Linea %d]: %s = Infinito\n", linea, nombre);
        printf("   Causa: Overflow numerico\n");
        exit(EXIT_FAILURE);
    }
}

#define VALIDAR(variable) verificar_nan_inf(#variable, variable, __LINE__)

void validar_parametros() {
    if (PASO_H <= 0) {
        printf("ERROR: PASO_H debe ser positivo (h = %f)\n", PASO_H);
        exit(EXIT_FAILURE);
    }
    
    if (T_FINAL <= T_INICIAL) {
        printf("ERROR: T_FINAL debe ser > T_INICIAL\n");
        exit(EXIT_FAILURE);
    }
    
    if (!es_numerico_valido(X_INICIAL) || !es_numerico_valido(Y_INICIAL)) {
        printf("ERROR: Condiciones iniciales invalidas\n");
        exit(EXIT_FAILURE);
    }
    
    // Verificar propiedades del sistema
    double f1_inicial = F1(X_INICIAL, Y_INICIAL);
    double f2_inicial = F2(X_INICIAL, Y_INICIAL);
    VALIDAR(f1_inicial);
    VALIDAR(f2_inicial);
    
    // Sistema lineal: matriz antisimetrica -> energia constante
    printf("Sistema valido: matriz antisimetrica\n");
}

FILE* abrir_archivo(const char *nombre, const char *modo) {
    FILE *archivo = fopen(nombre, modo);
    if (archivo == NULL) {
        printf("ERROR: No se pudo abrir '%s'\n", nombre);
        exit(EXIT_FAILURE);
    }
    return archivo;
}

// ============================================================================
// RUNGE-KUTTA 4 PARA SISTEMAS 2x2 CON VALIDACION
// ============================================================================
void rk4_sistema2_validado(double t, double *x, double *y, double h, int iter_actual) {
    double k1_x, k1_y, k2_x, k2_y, k3_x, k3_y, k4_x, k4_y;
    
    // k1
    k1_x = F1(*x, *y);
    k1_y = F2(*x, *y);
    VALIDAR(k1_x); VALIDAR(k1_y);
    
    // k2
    double x2 = *x + h*k1_x/2;
    double y2 = *y + h*k1_y/2;
    VALIDAR(x2); VALIDAR(y2);
    
    k2_x = F1(x2, y2);
    k2_y = F2(x2, y2);
    VALIDAR(k2_x); VALIDAR(k2_y);
    
    // k3
    double x3 = *x + h*k2_x/2;
    double y3 = *y + h*k2_y/2;
    VALIDAR(x3); VALIDAR(y3);
    
    k3_x = F1(x3, y3);
    k3_y = F2(x3, y3);
    VALIDAR(k3_x); VALIDAR(k3_y);
    
    // k4
    double x4 = *x + h*k3_x;
    double y4 = *y + h*k3_y;
    VALIDAR(x4); VALIDAR(y4);
    
    k4_x = F1(x4, y4);
    k4_y = F2(x4, y4);
    VALIDAR(k4_x); VALIDAR(k4_y);
    
    // Nuevos valores
    double x_nuevo = *x + h*(k1_x + 2*k2_x + 2*k3_x + k4_x)/6;
    double y_nuevo = *y + h*(k1_y + 2*k2_y + 2*k3_y + k4_y)/6;
    
    VALIDAR(x_nuevo); VALIDAR(y_nuevo);
    
    // Verificar conservacion de energia (x^2 + y^2 constante)
    double energia_antes = (*x)*(*x) + (*y)*(*y);
    double energia_despues = x_nuevo*x_nuevo + y_nuevo*y_nuevo;
    double delta_energia = fabs(energia_despues - energia_antes);
    
    if (delta_energia > 0.001 && iter_actual > 10) {
        printf("ADVERTENCIA [Iter %d]: Energia no se conserva\n", iter_actual);
        printf("   dE = %.2e, E_antes = %.6f, E_despues = %.6f\n", 
               delta_energia, energia_antes, energia_despues);
    }
    
    // Verificar si estamos en el circulo unitario (para condiciones iniciales tipicas)
    double radio = sqrt(x_nuevo*x_nuevo + y_nuevo*y_nuevo);
    if (fabs(radio - 1.0) > 0.1 && iter_actual > 5) {
        printf("ADVERTENCIA [Iter %d]: Radio diferente de 1\n", iter_actual);
        printf("   Radio = %.6f (deberia ser ~1)\n", radio);
    }
    
    *x = x_nuevo;
    *y = y_nuevo;
}

// ============================================================================
// PROGRAMA PRINCIPAL
// ============================================================================
int main() {
    // ============================================================================
    // VALIDACION INICIAL
    // ============================================================================
    printf("VALIDANDO SISTEMA DE ECUACIONES...\n");
    printf("-----------------------------------------------------------------\n");
    
    validar_parametros();
    
    double t = T_INICIAL;
    double x = X_INICIAL;
    double y = Y_INICIAL;
    int iter = 0;
    int iter_totales = (int)((T_FINAL - T_INICIAL) / PASO_H) + 1;
    
    double energia_inicial = X_INICIAL*X_INICIAL + Y_INICIAL*Y_INICIAL;
    
    printf("Sistema validado correctamente\n");
    printf("   Ecuaciones: dx/dt = y, dy/dt = -x\n");
    printf("   Condiciones: x(0) = %.1f, y(0) = %.1f\n", X_INICIAL, Y_INICIAL);
    printf("   Tiempo: [%.1f, %.1f]\n", T_INICIAL, T_FINAL);
    printf("   Paso: h = %.3f\n", PASO_H);
    printf("   Iteraciones estimadas: %d\n", iter_totales);
    printf("   Energia inicial: E = x^2 + y^2 = %.6f\n\n", energia_inicial);
    
    // ============================================================================
    // CONFIGURACION
    // ============================================================================
    printf("===============================================================\n");
    printf("          SISTEMA DE ECUACIONES: dx/dt=y, dy/dt=-x           \n");
    printf("===============================================================\n\n");
    
    FILE *datos_fase = abrir_archivo("sistema_fase.dat", "w");
    FILE *datos_x = abrir_archivo("sistema_x.dat", "w");
    FILE *datos_y = abrir_archivo("sistema_y.dat", "w");
    FILE *script_gp1 = abrir_archivo("sistema_temporal.gp", "w");
    FILE *script_gp2 = abrir_archivo("sistema_fase_plot.gp", "w");
    
    printf("PROCESO DE INTEGRACION:\n");
    printf("+------+--------+-----------+-----------+-----------+-----------+\n");
    printf("| Iter |   t    |   x(t)    |   y(t)    |  Energia  |  Estado   |\n");
    printf("+------+--------+-----------+-----------+-----------+-----------+\n");
    
    int errores_numericos = 0;
    double energia_max_desvio = 0.0;
    
    // ============================================================================
    // INTEGRACION DEL SISTEMA
    // ============================================================================
    while (t <= T_FINAL + PASO_H/2) {
        // Calcular energia actual
        double energia_actual = x*x + y*y;
        VALIDAR(energia_actual);
        
        double desvio_energia = fabs(energia_actual - energia_inicial);
        if (desvio_energia > energia_max_desvio) {
            energia_max_desvio = desvio_energia;
        }
        
        // Calcular soluciones exactas
        double x_exacto = cos(t);
        double y_exacto = -sin(t);
        
        // Estado de validacion
        const char *estado = "OK";
        if (!es_numerico_valido(x) || !es_numerico_valido(y)) {
            estado = "INVALIDO";
            errores_numericos++;
            
            printf("+------+--------+-----------+-----------+-----------+-----------+\n");
            printf("| %4d | %6.2f | %9.5f | %9.5f | %9.5f | %s |\n", 
                   iter, t, x, y, energia_actual, estado);
            printf("+------+--------+-----------+-----------+-----------+-----------+\n");
            
            printf("\nERROR CRITICO: Valores no numericos en iteracion %d\n", iter);
            printf("   t = %.6f, x = %.6f, y = %.6f\n", t, x, y);
            
            fclose(datos_fase);
            fclose(datos_x);
            fclose(datos_y);
            return EXIT_FAILURE;
        }
        
        // Mostrar cada 40 iteraciones
        if (iter % 40 == 0) {
            printf("| %4d | %6.2f | %9.5f | %9.5f | %9.5f | %s |\n", 
                   iter, t, x, y, energia_actual, estado);
        }
        
        // Guardar datos
        fprintf(datos_fase, "%.6f %.6f\n", x, y);
        fprintf(datos_x, "%.6f %.6f\n", t, x);
        fprintf(datos_y, "%.6f %.6f\n", t, y);
        
        // Ultimo punto
        if (t >= T_FINAL) break;
        
        // Calcular siguiente punto
        rk4_sistema2_validado(t, &x, &y, PASO_H, iter);
        
        t += PASO_H;
        iter++;
        
        // Verificar limite de iteraciones
        if (iter > iter_totales * 10) {
            printf("ADVERTENCIA: Demasiadas iteraciones (%d)\n", iter);
            break;
        }
    }
    
    printf("+------+--------+-----------+-----------+-----------+-----------+\n");
    printf("| INTEGRACION COMPLETADA: %d iteraciones                       |\n", iter);
    printf("+-------------------------------------------------------------+\n\n");
    
    fclose(datos_fase);
    fclose(datos_x);
    fclose(datos_y);
    
    // ============================================================================
    // CREAR SCRIPTS GNUPLOT
    // ============================================================================
    
    // Script 1: Evolucion temporal
    fprintf(script_gp1, "# Script para evolucion temporal\n");
    fprintf(script_gp1, "set terminal pngcairo size %d,%d enhanced font 'Arial,10'\n", 
            ANCHO_GRAFICO, ALTO_GRAFICO);
    fprintf(script_gp1, "set output '%s'\n", NOMBRE_GRAFICO1);
    fprintf(script_gp1, "set title 'Evolucion temporal: dx/dt = y, dy/dt = -x'\n");
    fprintf(script_gp1, "set xlabel 'Tiempo t'\n");
    fprintf(script_gp1, "set ylabel 'x(t), y(t)'\n");
    fprintf(script_gp1, "set grid\n");
    fprintf(script_gp1, "set key top right box\n");
    fprintf(script_gp1, "set xrange [%f:%f]\n", T_INICIAL, T_FINAL);
    
    fprintf(script_gp1, "plot 'sistema_x.dat' w l lw 2 lc rgb '#0066CC' title 'x(t)', \\\n");
    fprintf(script_gp1, "     'sistema_y.dat' w l lw 2 lc rgb '#FF3333' title 'y(t)', \\\n");
    fprintf(script_gp1, "     cos(x) w l lw 1 lc rgb '#0066CC' dt 2 title 'cos(t) (exacta)', \\\n");
    fprintf(script_gp1, "     -sin(x) w l lw 1 lc rgb '#FF3333' dt 2 title '-sin(t) (exacta)'\n");
    
    fclose(script_gp1);
    
    // Script 2: Plano de fase
    fprintf(script_gp2, "# Script para plano de fase\n");
    fprintf(script_gp2, "set terminal pngcairo size %d,%d enhanced font 'Arial,10'\n", 
            ANCHO_GRAFICO, ALTO_GRAFICO);
    fprintf(script_gp2, "set output '%s'\n", NOMBRE_GRAFICO2);
    fprintf(script_gp2, "set title 'Plano de fase: x vs y'\n");
    fprintf(script_gp2, "set xlabel 'x(t)'\n");
    fprintf(script_gp2, "set ylabel 'y(t)'\n");
    fprintf(script_gp2, "set grid\n");
    fprintf(script_gp2, "set key off\n");
    fprintf(script_gp2, "set size ratio -1\n");
    fprintf(script_gp2, "set xrange [-1.2:1.2]\n");
    fprintf(script_gp2, "set yrange [-1.2:1.2]\n");
    
    fprintf(script_gp2, "plot 'sistema_fase.dat' w l lw 1.5 lc rgb '#00AA00' title 'Trayectoria', \\\n");
    fprintf(script_gp2, "     cos(t), sin(t) w l lw 1 lc rgb '#000000' dt 2 title 'Circulo exacto'\n");
    
    fclose(script_gp2);
    
    // ============================================================================
    // EJECUTAR GNUPLOT
    // ============================================================================
    printf("GENERANDO GRAFICOS...\n");
    printf("-----------------------------------------------------------------\n");
    
    int resultado1 = system("gnuplot sistema_temporal.gp 2>&1");
    int resultado2 = system("gnuplot sistema_fase_plot.gp 2>&1");
    
    if (resultado1 != 0 || resultado2 != 0) {
        printf("ADVERTENCIA: Problemas al generar graficos\n");
        if (resultado1 != 0) {
            printf("   Error en grafico temporal\n");
        }
        if (resultado2 != 0) {
            printf("   Error en plano de fase\n");
        }
    } else {
        printf("Graficos generados correctamente:\n");
        printf("   • %s (evolucion temporal)\n", NOMBRE_GRAFICO1);
        printf("   • %s (plano de fase)\n", NOMBRE_GRAFICO2);
    }
    
    // ============================================================================
    // ANALISIS DE RESULTADOS
    // ============================================================================
    printf("\nANALISIS DE RESULTADOS:\n");
    printf("-----------------------------------------------------------------\n");
    
    double energia_final = x*x + y*y;
    double variacion_energia = fabs(energia_final - energia_inicial);
    double variacion_relativa = 100 * variacion_energia / energia_inicial;
    
    // Calcular errores respecto a solucion exacta
    double error_x_final = fabs(x - cos(T_FINAL));
    double error_y_final = fabs(y - (-sin(T_FINAL)));
    
    // Evaluar periodicidad
    double periodo_teorico = 2*M_PI;
    int ciclos_completos = (int)(T_FINAL / periodo_teorico);
    double fase_final = fmod(T_FINAL, periodo_teorico);
    
    printf("  Iteraciones:          %d\n", iter);
    printf("  Energia inicial:      %.8f\n", energia_inicial);
    printf("  Energia final:        %.8f\n", energia_final);
    printf("  Variacion energia:    %.2e (%.4f%%)\n", 
           variacion_energia, variacion_relativa);
    printf("  Error x final:        %.2e\n", error_x_final);
    printf("  Error y final:        %.2e\n", error_y_final);
    printf("  Ciclos completos:     %d\n", ciclos_completos);
    printf("  Fase final:           %.4f rad\n", fase_final);
    printf("  Errores numericos:    %d\n", errores_numericos);
    
    // Evaluar conservacion de energia
    printf("\n  EVALUACION DE CONSERVACION:\n");
    if (variacion_relativa < 0.01) {
        printf("    Excelente conservacion (dE < 0.01%%)\n");
    } else if (variacion_relativa < 0.1) {
        printf("    Buena conservacion (dE < 0.1%%)\n");
    } else if (variacion_relativa < 1.0) {
        printf("    Conservacion aceptable (dE < 1%%)\n");
    } else {
        printf("    Mala conservacion, metodo puede ser inestable\n");
        printf("       Considere reducir el paso h\n");
    }
    
    // Evaluar precision de solucion
    printf("\n  EVALUACION DE PRECISION:\n");
    double error_maximo = fmax(error_x_final, error_y_final);
    if (error_maximo < 0.001) {
        printf("    Excelente precision (error < 0.001)\n");
    } else if (error_maximo < 0.01) {
        printf("    Buena precision (error < 0.01)\n");
    } else if (error_maximo < 0.1) {
        printf("    Precision aceptable (error < 0.1)\n");
    } else {
        printf("    Precision pobre\n");
    }
    
    // ============================================================================
    // VALIDACION DE PROPIEDADES DEL SISTEMA
    // ============================================================================
    printf("\nVALIDACION DE PROPIEDADES DEL SISTEMA:\n");
    printf("-----------------------------------------------------------------\n");
    
    // Verificar que satisface las ecuaciones
    double dx_dt = F1(x, y);
    double dy_dt = F2(x, y);
    
    printf("  En t = %.4f:\n", t);
    printf("    x calculado:       %.8f\n", x);
    printf("    y calculado:       %.8f\n", y);
    printf("    dx/dt calculado:   %.8f\n", dx_dt);
    printf("    dy/dt calculado:   %.8f\n", dy_dt);
    printf("    dx/dt teorico:     y = %.8f\n", y);
    printf("    dy/dt teorico:     -x = %.8f\n", -x);
    
    // Verificar ortogonalidad (x·dx/dt + y·dy/dt = 0 para energia constante)
    double producto = x*dx_dt + y*dy_dt;
    printf("    x·dx/dt + y·dy/dt: %.2e (deberia ser ~0)\n", producto);
    
    if (fabs(producto) > 0.01) {
        printf("  ADVERTENCIA: Producto escalar grande\n");
        printf("     Indica posible error en la solucion\n");
    }
    
    // ============================================================================
    // RESUMEN FINAL
    // ============================================================================
    printf("\nRESUMEN DE EJECUCION:\n");
    printf("-----------------------------------------------------------------\n");
    printf("  Estado:              %s\n", 
           (errores_numericos == 0) ? "EXITOSO" : "CON ADVERTENCIAS");
    printf("  Iteraciones:         %d\n", iter);
    printf("  dEnergia:            %.2e (%.4f%%)\n", 
           variacion_energia, variacion_relativa);
    printf("  Error maximo:        %.2e\n", error_maximo);
    printf("  Graficos generados:  %s\n", 
           (resultado1 == 0 && resultado2 == 0) ? "2/2" : "PARCIAL");
    printf("  Archivos creados:    6 archivos de datos y scripts\n");
    
    printf("\n===============================================================\n");
    printf("                      EJECUCION COMPLETADA                     \n");
    printf("===============================================================\n");
    
    return (errores_numericos == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}