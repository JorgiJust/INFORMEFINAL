// 1_newton_raphson.c
// Metodo de Newton-Raphson con validaciones robustas

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <errno.h>

// ============================================================================

// ============================================================================
#define FUNCION(x)          ((x)*(x)*(x) - 2*(x) - 5)
#define DERIVADA(x)         (3*(x)*(x) - 2)
#define X_INICIAL           2.0
#define TOLERANCIA          1e-6
#define MAX_ITER            100
#define GRAFICO_INICIO      -3.0
#define GRAFICO_FIN         5.0
#define GRAFICO_PASO        0.1
#define NOMBRE_GRAFICO      "newton_grafico.png"
#define ANCHO_GRAFICO       800
#define ALTO_GRAFICO        600
// ============================================================================

// ============================================================================

// ============================================================================
int es_numerico_valido(double valor) {
    return !(isnan(valor) || isinf(valor) || fabs(valor) > 1e100);
}

void verificar_nan_inf(const char *nombre, double valor, int linea) {
    if (isnan(valor)) {
        printf(" ERROR en linea %d: %s = NaN (Not a Number)\n", linea, nombre);
        exit(EXIT_FAILURE);
    }
    if (isinf(valor)) {
        printf(" ERROR en linea %d: %s = Infinito\n", linea, nombre);
        exit(EXIT_FAILURE);
    }
    if (!es_numerico_valido(valor)) {
        printf(" ERROR en linea %d: %s = Valor numerico invalido\n", linea, nombre);
        exit(EXIT_FAILURE);
    }
}

#define VALIDAR(variable) verificar_nan_inf(#variable, variable, __LINE__)

FILE* abrir_archivo(const char *nombre, const char *modo) {
    FILE *archivo = fopen(nombre, modo);
    if (archivo == NULL) {
        printf(" ERROR: No se pudo abrir archivo '%s'\n", nombre);
        printf("   Verifique permisos o espacio en disco\n");
        exit(EXIT_FAILURE);
    }
    return archivo;
}

// ============================================================================
// FUNCIONES PRINCIPALES
// ============================================================================
void generar_datos_funcion() {
    FILE *func = abrir_archivo("funcion.dat", "w");
    fprintf(func, "# x f(x)\n");
    
    for (double xi = GRAFICO_INICIO; xi <= GRAFICO_FIN; xi += GRAFICO_PASO) {
        double fx = FUNCION(xi);
        VALIDAR(fx);
        fprintf(func, "%.3f %.3f\n", xi, fx);
    }
    fclose(func);
}

void crear_script_gnuplot(double raiz) {
    FILE *gp = abrir_archivo("newton_plot.gp", "w");
    
    fprintf(gp, "# Script Gnuplot para Newton-Raphson\n");
    fprintf(gp, "set terminal pngcairo size %d,%d enhanced font 'Arial,10'\n", 
            ANCHO_GRAFICO, ALTO_GRAFICO);
    fprintf(gp, "set output '%s'\n", NOMBRE_GRAFICO);
    fprintf(gp, "set title 'Metodo de Newton-Raphson: f(x) = x^3 - 2x - 5'\n");
    fprintf(gp, "set xlabel 'x'\n");
    fprintf(gp, "set ylabel 'f(x)'\n");
    fprintf(gp, "set grid\n");
    fprintf(gp, "set key top left box\n");
    fprintf(gp, "set zeroaxis lt -1\n\n");
    
    fprintf(gp, "plot 'funcion.dat' with lines lw 2 lc rgb 'blue' title 'f(x)', \\\n");
    fprintf(gp, "     0 with lines lc rgb 'black' notitle, \\\n");
    fprintf(gp, "     'iteraciones.dat' using 2:3 with points \\\n");
    fprintf(gp, "        pt 7 ps 1.5 lc rgb 'red' title 'Iteraciones', \\\n");
    fprintf(gp, "     %lf, 0 with points pt 9 ps 2 lc rgb 'green' title 'Raiz: %.6f'\n", 
            raiz, raiz);
    
    fclose(gp);
}

int ejecutar_gnuplot() {
    printf("\n Generando grafico...\n");
    int resultado = system("gnuplot newton_plot.gp 2>&1");
    
    if (resultado != 0) {
        printf(" ADVERTENCIA: Gnuplot encontro problemas\n");
        printf("   Verifique que Gnuplot este instalado: gnuplot --version\n");
        printf("   Puede generar el grafico manualmente con:\n");
        printf("   gnuplot newton_plot.gp\n");
        return 0;
    }
    
    printf(" Grafico generado exitosamente: %s\n", NOMBRE_GRAFICO);
    return 1;
}

int main() {
    double x = X_INICIAL, x_nuevo, error;
    int iter = 0;
    
    // ============================================================================
    // VALIDACION INICIAL DE PARAMS
    // ============================================================================
    printf(" Validando parametros iniciales...\n");
    
    if (!es_numerico_valido(X_INICIAL)) {
        printf("ERROR: Valor inicial X_INICIAL invalido: %f\n", X_INICIAL);
        return EXIT_FAILURE;
    }
    
    if (TOLERANCIA <= 0) {
        printf("ERROR: TOLERANCIA debe ser positiva: %e\n", TOLERANCIA);
        return EXIT_FAILURE;
    }
    
    if (MAX_ITER <= 0) {
        printf("ERROR: MAX_ITER debe ser positivo: %d\n", MAX_ITER);
        return EXIT_FAILURE;
    }
    
    double fx_inicial = FUNCION(x);
    double dfx_inicial = DERIVADA(x);
    
    VALIDAR(fx_inicial);
    VALIDAR(dfx_inicial);
    
    printf("Parametros validados correctamente\n\n");
    
    // ============================================================================

    // ============================================================================
    printf(" METODO DE NEWTON-RAPHSON \n\n");
    
    printf("CONFIGURACION:\n");
    printf("  Funcion:          f(x) = x3 - 2x - 5\n");
    printf("  Valor inicial:    x0 = %.1f, f(x0) = %.3f\n", X_INICIAL, fx_inicial);
    printf("  Derivada inicial: f'(x0) = %.3f\n", dfx_inicial);
    printf("  Tolerancia:       %.1e\n", TOLERANCIA);
    printf("  Max iteraciones:  %d\n\n", MAX_ITER);
    
    // Archivos
    FILE *datos = abrir_archivo("iteraciones.dat", "w");
    fprintf(datos, "# iter x f(x) error\n");
    
    generar_datos_funcion();
    
    printf("PROCESO DE CALCULO:\n");
    printf("+-----+-----------+-----------+-----------+-----------+\n");
    printf("| Iter|    x      |   f(x)    |  f'(x)    |  Error    |\n");
    printf("+-----+-----------+-----------+-----------+-----------+\n");
    
    // ============================================================================
    // NEWTON-RAPHSON CON VALIDACIONES
    // ============================================================================
    do {
        double fx = FUNCION(x);
        double dfx = DERIVADA(x);
        
        VALIDAR(fx);
        VALIDAR(dfx);
        
        // Validacion de derivada
        if (fabs(dfx) < 1e-15) {
            printf("+-----+-----------+-----------+-----------+-----------+\n");
            printf("| ERROR CRITICO: Derivada cero (%.2e)                |\n", dfx);
            printf("|   en x = %.6f                                        |\n", x);
            printf("|   f(x) = %.6f                                        |\n", fx);
            printf("|   El metodo no puede continuar                        |\n");
            printf("+-----------------------------------------------------+\n");
            fclose(datos);
            return EXIT_FAILURE;
        }
        
        // nueva aproximacion
        x_nuevo = x - fx / dfx;
        VALIDAR(x_nuevo);
        
        error = fabs(x_nuevo - x);
        VALIDAR(error);
        
        //  divergencia
        if (error > 1e10 && iter > 5) {
            printf("+-----+-----------+-----------+-----------+-----------+\n");
            printf("| ADVERTENCIA: Posible divergencia                   |\n");
            printf("|   Error creciente: %.2e                            |\n", error);
            printf("|   Considere cambiar el valor inicial                 |\n");
            printf("+-----------------------------------------------------+\n");
            break;
        }
        
        // Mostrar y guardar
        printf("| %3d | %9.6f | %9.6f | %9.6f | %9.6f |\n", 
               iter, x, fx, dfx, error);
        
        fprintf(datos, "%d %.6f %.6f %.6f\n", iter, x, fx, error);
        
        // Actualizar
        x = x_nuevo;
        iter++;
        
        // Verificar convergencia
        if (error < TOLERANCIA) {
            printf("+-----+-----------+-----------+-----------+-----------+\n");
            printf("| CONVERGENCIA ALCANZADA                            |\n");
            printf("|   Error final: %.2e < Tolerancia: %.2e           |\n", 
                   error, TOLERANCIA);
            printf("+-----------------------------------------------------+\n\n");
            break;
        }
        
        // Verificar max de iteraciones
        if (iter >= MAX_ITER) {
            printf("+-----+-----------+-----------+-----------+-----------+\n");
            printf("| LIMITE DE ITERACIONES ALCANZADO                    |\n");
            printf("|   No se alcanzo la tolerancia en %d iteraciones     |\n", MAX_ITER);
            printf("|   Ultimo error: %.2e                            |\n", error);
            printf("+-----------------------------------------------------+\n\n");
            break;
        }
        
    } while (1);
    
    fclose(datos);
    
    // Validar resultado final
    double fx_final = FUNCION(x);
    VALIDAR(fx_final);
    
    if (fabs(fx_final) > 0.1) {
        printf(" ADVERTENCIA: Valor de funcion en raiz es alto: %.2e\n", fx_final);
        printf("   La raiz podria no ser precisa\n");
    }
    
    // ============================================================================
    // GENERAR 
    // ============================================================================
    crear_script_gnuplot(x);
    int grafico_ok = ejecutar_gnuplot();
    
    // ============================================================================
    // RESULTADOS FINALES
    // ============================================================================
    printf("\n RESULTADOS FINALES:\n");
    printf("-------------------------------------------------------------\n");
    printf("  Raiz aproximada:  x = %.8f\n", x);
    printf("  f(raiz) =         %.2e\n", fx_final);
    printf("  Iteraciones:      %d de %d\n", iter, MAX_ITER);
    printf("  Error final:      %.2e (Tolerancia: %.2e)\n", error, TOLERANCIA);
    printf("  Estado:           %s\n", 
           (error < TOLERANCIA) ? "CONVERGENCIA" : "ITERACIONES MAXIMAS");
    printf("  Grafico:          %s\n", 
           grafico_ok ? "GENERADO CORRECTAMENTE" : "NO SE PUDO GENERAR");
    
    printf("\n ARCHIVOS GENERADOS:\n");
    printf("-------------------------------------------------------------\n");
    printf("  - iteraciones.dat   -> %d iteraciones guardadas\n", iter);
    printf("  - funcion.dat       -> Puntos para graficar\n");
    printf("  - newton_plot.gp    -> Script de Gnuplot\n");
    if (grafico_ok) {
        printf("  - %s -> Grafico final\n", NOMBRE_GRAFICO);
    }
    
    return EXIT_SUCCESS;
}