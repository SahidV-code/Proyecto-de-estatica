/**
 * ============================================================
 *  ANALISIS DE ARMADURAS PLANAS - METODO DE NODOS
 *  Practica ICC - CETYS Universidad
 *  Docente: Ulises Alejandro Del Moral Phelps
 * ============================================================
 *
 *  Compilar:  g++ -std=c++17 -O2 -o armaduras armaduras.cpp
 *  Ejecutar:  ./armaduras
 *
 *  El usuario puede:
 *    1) Ingresar una armadura manualmente (nodos, barras, apoyos, cargas)
 *    2) Cargar un caso de estudio predefinido (AC.1 - AC.4)
 *    3) Ejecutar todos los casos predefinidos de una vez
 *    4) Salir
 *
 *  Convenios de signo:
 *    - Fuerza axial positiva = TENSION
 *    - Fuerza axial negativa = COMPRESION
 *    - Cargas: Fx positivo = derecha, Fy positivo = arriba
 *
 *  Apoyos soportados:
 *    - PASADOR  (pin)    : reacciones Rx y Ry
 *    - RODILLO_Y         : reaccion Ry  (rueda sobre superficie horizontal)
 *    - RODILLO_X         : reaccion Rx  (rueda sobre superficie vertical)
 * ============================================================
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <stdexcept>
#include <map>
#include <limits>

 // ─────────────────────────────────────────────────────────────
 //  Constantes
 // ─────────────────────────────────────────────────────────────
static const double ZERO_TOL = 1e-6;
static const double SING_TOL = 1e-12;

// ─────────────────────────────────────────────────────────────
//  Estructuras de datos
// ─────────────────────────────────────────────────────────────

struct Nodo {
    std::string nombre;
    double x, y;
    double fx = 0.0;
    double fy = 0.0;
};

struct Barra {
    std::string nombre;
    int i, j;
};

enum class TipoApoyo { PASADOR, RODILLO_X, RODILLO_Y };

struct Apoyo {
    int nodo_idx;
    TipoApoyo tipo;
};

struct ResultadoBarra {
    std::string nombre;
    double fuerza;
    std::string estado;
};

struct ResultadoArmadura {
    std::vector<ResultadoBarra> barras;
    std::map<std::string, double> reacciones;
    std::vector<double> residual_nodo;
    double residual_max;
    bool isostatica;
    int m, r, j;
};

// ─────────────────────────────────────────────────────────────
//  Algebra lineal: eliminacion gaussiana con pivoteo parcial
// ─────────────────────────────────────────────────────────────

std::vector<double> resolver_sistema(std::vector<std::vector<double>> A,
    std::vector<double> b)
{
    int n = (int)b.size();
    for (int col = 0; col < n; ++col) {
        int pivot_fila = col;
        double max_val = std::abs(A[col][col]);
        for (int fila = col + 1; fila < n; ++fila)
            if (std::abs(A[fila][col]) > max_val) {
                max_val = std::abs(A[fila][col]);
                pivot_fila = fila;
            }
        if (max_val < SING_TOL)
            throw std::runtime_error(
                "Sistema singular. Verifique geometria y apoyos.");
        std::swap(A[col], A[pivot_fila]);
        std::swap(b[col], b[pivot_fila]);
        for (int fila = col + 1; fila < n; ++fila) {
            double factor = A[fila][col] / A[col][col];
            for (int k = col; k < n; ++k)
                A[fila][k] -= factor * A[col][k];
            b[fila] -= factor * b[col];
        }
    }
    std::vector<double> x(n, 0.0);
    for (int i = n - 1; i >= 0; --i) {
        x[i] = b[i];
        for (int k = i + 1; k < n; ++k) x[i] -= A[i][k] * x[k];
        x[i] /= A[i][i];
    }
    return x;
}

// ─────────────────────────────────────────────────────────────
//  Motor principal
// ─────────────────────────────────────────────────────────────

ResultadoArmadura resolver_armadura(const std::vector<Nodo>& nodos,
    const std::vector<Barra>& barras,
    const std::vector<Apoyo>& apoyos)
{
    ResultadoArmadura res;
    int j = (int)nodos.size();
    int m = (int)barras.size();

    struct Reaccion { int nodo_idx; char dir; std::string nombre; };
    std::vector<Reaccion> reacciones_list;
    for (const auto& ap : apoyos) {
        const std::string& nn = nodos[ap.nodo_idx].nombre;
        if (ap.tipo == TipoApoyo::PASADOR) {
            reacciones_list.push_back({ ap.nodo_idx, 'x', "Rx_" + nn });
            reacciones_list.push_back({ ap.nodo_idx, 'y', "Ry_" + nn });
        }
        else if (ap.tipo == TipoApoyo::RODILLO_Y) {
            reacciones_list.push_back({ ap.nodo_idx, 'y', "Ry_" + nn });
        }
        else {
            reacciones_list.push_back({ ap.nodo_idx, 'x', "Rx_" + nn });
        }
    }

    int r = (int)reacciones_list.size();
    res.m = m; res.r = r; res.j = j;
    res.isostatica = (m + r == 2 * j);

    std::cout << "\n  Verificacion isostatica: m=" << m
        << "  r=" << r << "  j=" << j
        << "  ->  m+r=" << (m + r) << "  2j=" << (2 * j);
    std::cout << (res.isostatica ? "  [OK]\n" : "  [!! NO CUMPLE]\n");

    if (m + r != 2 * j)
        throw std::runtime_error("m+r != 2j. El sistema no es cuadrado.");

    int neq = 2 * j;
    int nunk = m + r;
    std::vector<std::vector<double>> A(neq, std::vector<double>(nunk, 0.0));
    std::vector<double> b(neq, 0.0);

    for (int k = 0; k < m; ++k) {
        int ni = barras[k].i, nj = barras[k].j;
        double dx = nodos[nj].x - nodos[ni].x;
        double dy = nodos[nj].y - nodos[ni].y;
        double L = std::sqrt(dx * dx + dy * dy);
        if (L < ZERO_TOL)
            throw std::runtime_error("Barra '" + barras[k].nombre + "' tiene longitud cero.");
        double c = dx / L, s = dy / L;
        A[2 * ni][k] += c;  A[2 * ni + 1][k] += s;
        A[2 * nj][k] -= c;  A[2 * nj + 1][k] -= s;
    }
    for (int ri = 0; ri < r; ++ri) {
        int ni = reacciones_list[ri].nodo_idx;
        int col = m + ri;
        if (reacciones_list[ri].dir == 'x') A[2 * ni][col] += 1.0;
        else                                 A[2 * ni + 1][col] += 1.0;
    }
    for (int i = 0; i < j; ++i) {
        b[2 * i] = -nodos[i].fx;
        b[2 * i + 1] = -nodos[i].fy;
    }

    std::vector<double> x = resolver_sistema(A, b);

    for (int k = 0; k < m; ++k) {
        ResultadoBarra rb;
        rb.nombre = barras[k].nombre;
        rb.fuerza = x[k];
        if (x[k] > ZERO_TOL) rb.estado = "TENSION";
        else if (x[k] < -ZERO_TOL) rb.estado = "COMPRESION";
        else                        rb.estado = "FUERZA CERO";
        res.barras.push_back(rb);
    }
    for (int ri = 0; ri < r; ++ri)
        res.reacciones[reacciones_list[ri].nombre] = x[m + ri];

    res.residual_max = 0.0;
    for (int eq = 0; eq < neq; ++eq) {
        double val = 0.0;
        for (int col = 0; col < nunk; ++col) val += A[eq][col] * x[col];
        val -= b[eq];
        res.residual_nodo.push_back(std::abs(val));
        res.residual_max = std::max(res.residual_max, std::abs(val));
    }
    return res;
}

// ─────────────────────────────────────────────────────────────
//  Impresion de resultados
// ─────────────────────────────────────────────────────────────

void imprimir_resultado(const std::string& titulo,
    const ResultadoArmadura& res,
    const std::vector<Nodo>& nodos)
{
    const int W = 62;
    std::cout << "\n" << std::string(W, '=') << "\n";
    std::cout << "  " << titulo << "\n";
    std::cout << std::string(W, '=') << "\n";

    std::cout << "\n  FUERZAS INTERNAS EN BARRAS\n";
    std::cout << "  " << std::string(W - 2, '-') << "\n";
    std::cout << std::left
        << "  " << std::setw(10) << "Barra"
        << std::setw(18) << "Fuerza (kN)"
        << std::setw(14) << "Magnitud"
        << "Estado\n";
    std::cout << "  " << std::string(W - 2, '-') << "\n";
    for (const auto& b : res.barras)
        std::cout << std::fixed << std::setprecision(4)
        << "  " << std::left << std::setw(10) << b.nombre
        << std::right << std::setw(14) << b.fuerza
        << std::setw(14) << std::abs(b.fuerza)
        << "  " << b.estado << "\n";

    std::cout << "\n  REACCIONES DE APOYO\n";
    std::cout << "  " << std::string(W - 2, '-') << "\n";
    for (const auto& [nombre, val] : res.reacciones)
        std::cout << std::fixed << std::setprecision(4)
        << "  " << std::left << std::setw(14) << nombre
        << std::right << std::setw(14) << val << " kN\n";

    std::cout << "\n  RESIDUALES DE EQUILIBRIO POR NODO\n";
    std::cout << "  " << std::string(W - 2, '-') << "\n";
    for (int i = 0; i < (int)nodos.size(); ++i)
        std::cout << std::scientific << std::setprecision(3)
        << "  Nodo " << std::left << std::setw(6) << nodos[i].nombre
        << "  |res_x|=" << res.residual_nodo[2 * i]
        << "  |res_y|=" << res.residual_nodo[2 * i + 1] << "\n";
    std::cout << std::scientific << std::setprecision(3)
        << "\n  Residual maximo global: " << res.residual_max << "\n";
    std::cout << std::string(W, '=') << "\n";
}

// ─────────────────────────────────────────────────────────────
//  Exportar CSV
// ─────────────────────────────────────────────────────────────

void exportar_csv(const std::string& archivo, const std::string& titulo,
    const ResultadoArmadura& res, const std::vector<Nodo>& nodos)
{
    std::ofstream f(archivo);
    if (!f.is_open()) { std::cerr << "No se pudo crear: " << archivo << "\n"; return; }
    f << "CASO," << titulo << "\n";
    f << "m+r=2j," << (res.m + res.r) << "=" << (2 * res.j)
        << "," << (res.isostatica ? "CUMPLE" : "NO CUMPLE") << "\n\n";
    f << "FUERZAS INTERNAS\nBarra,Fuerza (kN),Magnitud (kN),Estado\n";
    for (const auto& b : res.barras)
        f << b.nombre << "," << std::fixed << std::setprecision(6)
        << b.fuerza << "," << std::abs(b.fuerza) << "," << b.estado << "\n";
    f << "\nREACCIONES\nReaccion,Valor (kN)\n";
    for (const auto& [n, v] : res.reacciones)
        f << n << "," << std::fixed << std::setprecision(6) << v << "\n";
    f << "\nRESIDUALES\nNodo,|res_x|,|res_y|\n";
    for (int i = 0; i < (int)nodos.size(); ++i)
        f << nodos[i].nombre << "," << std::scientific << std::setprecision(4)
        << res.residual_nodo[2 * i] << "," << res.residual_nodo[2 * i + 1] << "\n";
    f << "\nResidual maximo," << std::scientific << res.residual_max << "\n";
    f.close();
    std::cout << "  -> CSV exportado: " << archivo << "\n";
}

// ─────────────────────────────────────────────────────────────
//  Helpers
// ─────────────────────────────────────────────────────────────

int idx_nodo(const std::vector<Nodo>& nodos, const std::string& nombre) {
    for (int i = 0; i < (int)nodos.size(); ++i)
        if (nodos[i].nombre == nombre) return i;
    return -1;
}

Barra hacer_barra(const std::vector<Nodo>& nodos,
    const std::string& ni, const std::string& nj)
{
    Barra b;
    b.nombre = ni + "-" + nj;
    b.i = idx_nodo(nodos, ni);
    b.j = idx_nodo(nodos, nj);
    if (b.i < 0 || b.j < 0)
        throw std::runtime_error("Nodo no encontrado: " + ni + " o " + nj);
    return b;
}

const ResultadoBarra& barra_max_abs(const std::vector<ResultadoBarra>& barras) {
    return *std::max_element(barras.begin(), barras.end(),
        [](const ResultadoBarra& a, const ResultadoBarra& b) {
            return std::abs(a.fuerza) < std::abs(b.fuerza); });
}

// ─────────────────────────────────────────────────────────────
//  Leer string sin espacios
// ─────────────────────────────────────────────────────────────

std::string leer_string(const std::string& prompt) {
    std::string s;
    std::cout << prompt;
    std::cin >> s;
    return s;
}

double leer_double(const std::string& prompt) {
    double v;
    while (true) {
        std::cout << prompt;
        if (std::cin >> v) return v;
        std::cin.clear();
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        std::cout << "  [Error] Ingrese un numero valido.\n";
    }
}

int leer_int(const std::string& prompt, int min_val, int max_val) {
    int v;
    while (true) {
        std::cout << prompt;
        if (std::cin >> v && v >= min_val && v <= max_val) return v;
        std::cin.clear();
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        std::cout << "  [Error] Ingrese un entero entre " << min_val
            << " y " << max_val << ".\n";
    }
}

// ─────────────────────────────────────────────────────────────
//  ENTRADA MANUAL DE ARMADURA
// ─────────────────────────────────────────────────────────────

void ingresar_armadura_manual() {
    std::cout << "\n============================================================\n";
    std::cout << "  INGRESO MANUAL DE ARMADURA\n";
    std::cout << "============================================================\n";

    // ── Nodos ─────────────────────────────────────────────────
    int num_nodos = leer_int("\nNumero de nodos: ", 2, 100);
    std::vector<Nodo> nodos;
    std::cout << "\n  Ingrese coordenadas de cada nodo:\n";
    for (int i = 0; i < num_nodos; ++i) {
        Nodo nd;
        nd.nombre = leer_string("  Nombre del nodo " + std::to_string(i + 1) + ": ");
        nd.x = leer_double("    x (m): ");
        nd.y = leer_double("    y (m): ");
        nd.fx = 0.0; nd.fy = 0.0;
        nodos.push_back(nd);
    }

    // ── Barras ────────────────────────────────────────────────
    int num_barras = leer_int("\nNumero de barras: ", 1, 200);
    std::vector<Barra> barras;
    std::cout << "\n  Ingrese las barras (por nombre de nodo):\n";
    for (int k = 0; k < num_barras; ++k) {
        std::cout << "  Barra " << (k + 1) << ":\n";
        std::string ni = leer_string("    Nodo inicial: ");
        std::string nj = leer_string("    Nodo final  : ");
        try {
            barras.push_back(hacer_barra(nodos, ni, nj));
        }
        catch (const std::exception& e) {
            std::cout << "  [Error] " << e.what() << " - barra omitida.\n";
            --k;
        }
    }

    // ── Apoyos ────────────────────────────────────────────────
    int num_apoyos = leer_int("\nNumero de apoyos: ", 1, 20);
    std::vector<Apoyo> apoyos;
    std::cout << "\n  Tipos de apoyo:\n"
        << "    1) Pasador  (Rx + Ry)\n"
        << "    2) Rodillo horizontal - Ry  (rueda sobre piso)\n"
        << "    3) Rodillo vertical   - Rx  (rueda contra pared)\n\n";
    for (int a = 0; a < num_apoyos; ++a) {
        std::cout << "  Apoyo " << (a + 1) << ":\n";
        std::string nn = leer_string("    Nombre del nodo: ");
        int idx = idx_nodo(nodos, nn);
        if (idx < 0) {
            std::cout << "  [Error] Nodo no encontrado. Intente de nuevo.\n";
            --a; continue;
        }
        int tipo = leer_int("    Tipo (1/2/3): ", 1, 3);
        TipoApoyo ta = (tipo == 1) ? TipoApoyo::PASADOR
            : (tipo == 2) ? TipoApoyo::RODILLO_Y
            : TipoApoyo::RODILLO_X;
        apoyos.push_back({ idx, ta });
    }

    // ── Cargas externas ───────────────────────────────────────
    int num_cargas = leer_int("\nNumero de nodos con carga aplicada: ", 0, 50);
    for (int c = 0; c < num_cargas; ++c) {
        std::cout << "  Carga " << (c + 1) << ":\n";
        std::string nn = leer_string("    Nombre del nodo: ");
        int idx = idx_nodo(nodos, nn);
        if (idx < 0) {
            std::cout << "  [Error] Nodo no encontrado. Intente de nuevo.\n";
            --c; continue;
        }
        std::cout << "    (positivo = derecha/arriba, negativo = izquierda/abajo)\n";
        nodos[idx].fx = leer_double("    Fx (kN): ");
        nodos[idx].fy = leer_double("    Fy (kN): ");
    }

    // ── Resolver e imprimir ───────────────────────────────────
    try {
        auto res = resolver_armadura(nodos, barras, apoyos);
        imprimir_resultado("Armadura ingresada por usuario", res, nodos);

        std::cout << "\n  Desea exportar resultados a CSV? (1=Si / 0=No): ";
        int op; std::cin >> op;
        if (op == 1) {
            std::string arch = leer_string("  Nombre del archivo (sin extension): ");
            exportar_csv(arch + ".csv", "Armadura usuario", res, nodos);
        }
    }
    catch (const std::exception& e) {
        std::cout << "\n  [ERROR] " << e.what() << "\n";
    }
}

// ─────────────────────────────────────────────────────────────
//  CASOS PREDEFINIDOS
// ─────────────────────────────────────────────────────────────

void caso_AC1() {
    std::cout << "\n============================================================\n";
    std::cout << "  CASO AC.1 - ARMADURA TRIANGULAR BASICA\n";
    std::cout << "============================================================\n";

    auto build = [](double P) {
        std::vector<Nodo> nodos = {
            {"A",0,0,0,0}, {"B",3,4,0,-P}, {"C",6,0,0,0}
        };
        std::vector<Barra> barras = {
            hacer_barra(nodos,"A","B"),
            hacer_barra(nodos,"B","C"),
            hacer_barra(nodos,"A","C")
        };
        std::vector<Apoyo> apoyos = {
            {idx_nodo(nodos,"A"), TipoApoyo::PASADOR},
            {idx_nodo(nodos,"C"), TipoApoyo::RODILLO_Y}
        };
        return std::make_tuple(nodos, barras, apoyos);
        };

    double P = leer_double("\n  Carga PB (kN, hacia abajo positivo): ");
    auto [nodos, barras, apoyos] = build(P);
    auto res = resolver_armadura(nodos, barras, apoyos);
    imprimir_resultado("AC.1 | PB = " + std::to_string(P) + " kN", res, nodos);
    exportar_csv("AC1_resultados.csv", "AC.1", res, nodos);

    // Estudio parametrico
    std::cout << "\n  Ejecutar estudio parametrico PB = 10,20,30,40,50 kN? (1=Si/0=No): ";
    int op; std::cin >> op;
    if (op == 1) {
        std::cout << "\n  ESTUDIO PARAMETRICO AC.1\n";
        std::cout << "  " << std::string(58, '-') << "\n";
        std::cout << std::left
            << "  " << std::setw(10) << "PB (kN)"
            << std::setw(18) << "Max|F| barra"
            << std::setw(16) << "Fuerza (kN)"
            << "Estado\n";
        std::cout << "  " << std::string(58, '-') << "\n";
        std::ofstream fcsv("AC1_parametrico.csv");
        fcsv << "PB (kN),Barra max,Fuerza (kN),|Fuerza| (kN),Estado\n";
        for (double Pv : {10.0, 20.0, 30.0, 40.0, 50.0}) {
            auto [n2, b2, a2] = build(Pv);
            auto r2 = resolver_armadura(n2, b2, a2);
            const auto& bm = barra_max_abs(r2.barras);
            std::cout << std::fixed << std::setprecision(2)
                << "  " << std::setw(10) << Pv
                << std::setw(18) << bm.nombre
                << std::setw(16) << bm.fuerza
                << bm.estado << "\n";
            fcsv << std::fixed << std::setprecision(6)
                << Pv << "," << bm.nombre << "," << bm.fuerza << ","
                << std::abs(bm.fuerza) << "," << bm.estado << "\n";
        }
        fcsv.close();
        std::cout << "  -> CSV parametrico exportado: AC1_parametrico.csv\n";
    }
}

void caso_AC2() {
    std::cout << "\n============================================================\n";
    std::cout << "  CASO AC.2 - ARMADURA WARREN DE DOS CLAROS\n";
    std::cout << "============================================================\n";

    double PB = leer_double("\n  Carga PB (kN, hacia abajo positivo): ");
    double PD = leer_double("  Carga PD (kN, hacia abajo positivo): ");

    std::vector<Nodo> nodos = {
        {"A",0,0,0,0}, {"B",2,3,0,-PB}, {"C",4,0,0,0},
        {"D",6,3,0,-PD}, {"E",8,0,0,0}
    };
    std::vector<Barra> barras = {
        hacer_barra(nodos,"A","B"), hacer_barra(nodos,"B","C"),
        hacer_barra(nodos,"A","C"), hacer_barra(nodos,"B","D"),
        hacer_barra(nodos,"C","D"), hacer_barra(nodos,"D","E"),
        hacer_barra(nodos,"C","E")
    };
    std::vector<Apoyo> apoyos = {
        {idx_nodo(nodos,"A"), TipoApoyo::PASADOR},
        {idx_nodo(nodos,"E"), TipoApoyo::RODILLO_Y}
    };

    auto res = resolver_armadura(nodos, barras, apoyos);
    imprimir_resultado("AC.2 | Warren PB=" + std::to_string(PB)
        + " PD=" + std::to_string(PD) + " kN", res, nodos);
    exportar_csv("AC2_resultados.csv", "AC.2 Warren", res, nodos);

    // Ordenar por magnitud
    std::vector<ResultadoBarra> ord = res.barras;
    std::sort(ord.begin(), ord.end(),
        [](const ResultadoBarra& a, const ResultadoBarra& b) {
            return std::abs(a.fuerza) > std::abs(b.fuerza); });
    std::cout << "\n  BARRAS ORDENADAS POR MAGNITUD (mayor -> menor)\n";
    std::cout << "  " << std::string(54, '-') << "\n";
    std::cout << std::left
        << "  " << std::setw(4) << "Rk"
        << std::setw(12) << "Barra"
        << std::setw(18) << "Fuerza (kN)"
        << "Estado\n";
    std::cout << "  " << std::string(54, '-') << "\n";
    int rk = 1;
    for (const auto& b : ord)
        std::cout << std::fixed << std::setprecision(4)
        << "  " << std::setw(4) << rk++
        << std::setw(12) << b.nombre
        << std::setw(18) << b.fuerza
        << b.estado << "\n";
}

void caso_AC3() {
    std::cout << "\n============================================================\n";
    std::cout << "  CASO AC.3 - CARGA OBLICUA RESULTANTE\n";
    std::cout << "============================================================\n";

    double Px = leer_double("\n  Carga Px en B (kN, derecha positivo): ");
    double Py = leer_double("  Carga Py en B (kN, abajo positivo): ");

    std::vector<Nodo> nodos = {
        {"A",0,0,0,0}, {"B",3,4,Px,-Py}, {"C",6,0,0,0}, {"D",3,0,0,0}
    };
    std::vector<Barra> barras = {
        hacer_barra(nodos,"A","B"), hacer_barra(nodos,"B","D"),
        hacer_barra(nodos,"B","C"), hacer_barra(nodos,"A","D"),
        hacer_barra(nodos,"D","C")
    };
    std::vector<Apoyo> apoyos = {
        {idx_nodo(nodos,"A"), TipoApoyo::PASADOR},
        {idx_nodo(nodos,"C"), TipoApoyo::RODILLO_Y}
    };

    auto res = resolver_armadura(nodos, barras, apoyos);
    imprimir_resultado("AC.3 | Carga oblicua Px=" + std::to_string(Px)
        + " Py=" + std::to_string(Py) + " kN", res, nodos);
    exportar_csv("AC3_resultados.csv", "AC.3 Oblicua", res, nodos);
}

void caso_AC4() {
    std::cout << "\n============================================================\n";
    std::cout << "  CASO AC.4 - ARMADURA TIPO PRATT SIMPLIFICADA\n";
    std::cout << "============================================================\n";

    double P = leer_double("\n  Carga P en B, D y F (kN, hacia abajo positivo): ");

    auto build = [](double Pv) {
        std::vector<Nodo> nodos = {
            {"A", 0,0,0,0}, {"C", 3,0,0,0}, {"E", 6,0,0,0},
            {"G", 9,0,0,0}, {"I",12,0,0,0},
            {"B", 3,3,0,-Pv}, {"D", 6,3,0,-Pv}, {"F", 9,3,0,-Pv}
        };
        std::vector<Barra> barras = {
            hacer_barra(nodos,"A","C"), hacer_barra(nodos,"C","E"),
            hacer_barra(nodos,"E","G"), hacer_barra(nodos,"G","I"),
            hacer_barra(nodos,"B","D"), hacer_barra(nodos,"D","F"),
            hacer_barra(nodos,"C","B"), hacer_barra(nodos,"E","D"),
            hacer_barra(nodos,"G","F"), hacer_barra(nodos,"A","B"),
            hacer_barra(nodos,"F","I"), hacer_barra(nodos,"B","E"),
            hacer_barra(nodos,"D","G")
        };
        std::vector<Apoyo> apoyos = {
            {idx_nodo(nodos,"A"), TipoApoyo::PASADOR},
            {idx_nodo(nodos,"I"), TipoApoyo::RODILLO_Y}
        };
        return std::make_tuple(nodos, barras, apoyos);
        };

    {
        auto [nodos, barras, apoyos] = build(P);
        auto res = resolver_armadura(nodos, barras, apoyos);
        imprimir_resultado("AC.4 | Pratt P=" + std::to_string(P) + " kN", res, nodos);
        exportar_csv("AC4_resultados.csv", "AC.4 Pratt", res, nodos);
    }

    std::cout << "\n  Ejecutar estudio parametrico P = 5,10,15,20,25 kN? (1=Si/0=No): ";
    int op; std::cin >> op;
    if (op == 1) {
        std::cout << "\n  ESTUDIO PARAMETRICO AC.4 - Pratt\n";
        std::cout << "  " << std::string(70, '-') << "\n";
        std::cout << std::left
            << "  " << std::setw(10) << "P (kN)"
            << std::setw(16) << "Max tension"
            << std::setw(14) << "F (kN)"
            << std::setw(18) << "Max compresion"
            << "F (kN)\n";
        std::cout << "  " << std::string(70, '-') << "\n";

        std::ofstream fcsv("AC4_parametrico.csv");
        fcsv << "P (kN),Barra max tension,F_tension (kN),"
            "Barra max compresion,F_compresion (kN)\n";

        for (double Pv : {5.0, 10.0, 15.0, 20.0, 25.0}) {
            auto [n2, b2, a2] = build(Pv);
            auto r2 = resolver_armadura(n2, b2, a2);
            const ResultadoBarra* mt = nullptr, * mc = nullptr;
            for (const auto& b : r2.barras) {
                if (b.fuerza > ZERO_TOL && (!mt || b.fuerza > mt->fuerza)) mt = &b;
                if (b.fuerza < -ZERO_TOL && (!mc || b.fuerza < mc->fuerza)) mc = &b;
            }
            double ft = mt ? mt->fuerza : 0.0;
            double fc = mc ? mc->fuerza : 0.0;
            std::string nt = mt ? mt->nombre : "N/A";
            std::string nc = mc ? mc->nombre : "N/A";
            std::cout << std::fixed << std::setprecision(2)
                << "  " << std::setw(10) << Pv
                << std::setw(16) << nt << std::setw(14) << ft
                << std::setw(18) << nc << fc << "\n";
            fcsv << std::fixed << std::setprecision(6)
                << Pv << "," << nt << "," << ft << ","
                << nc << "," << fc << "\n";
        }
        fcsv.close();
        std::cout << "  -> CSV parametrico exportado: AC4_parametrico.csv\n";
        std::cout << "\n  NOTA: La respuesta es LINEAL porque A no depende de P;\n"
            << "  solo el vector b escala con P => F_k(aP) = a*F_k(P).\n";
    }
}

void ejecutar_todos() {
    std::cout << "\n  Ejecutando todos los casos con valores predeterminados...\n";

    // AC1 P=20
    {
        std::vector<Nodo> nodos = { {"A",0,0,0,0},{"B",3,4,0,-20},{"C",6,0,0,0} };
        std::vector<Barra> barras = { hacer_barra(nodos,"A","B"),hacer_barra(nodos,"B","C"),hacer_barra(nodos,"A","C") };
        std::vector<Apoyo> apoyos = { {idx_nodo(nodos,"A"),TipoApoyo::PASADOR},{idx_nodo(nodos,"C"),TipoApoyo::RODILLO_Y} };
        auto res = resolver_armadura(nodos, barras, apoyos);
        imprimir_resultado("AC.1 | PB=20 kN", res, nodos);
        exportar_csv("AC1_resultados.csv", "AC.1 PB=20kN", res, nodos);
    }
    // AC2 PB=15 PD=25
    {
        std::vector<Nodo> nodos = { {"A",0,0,0,0},{"B",2,3,0,-15},{"C",4,0,0,0},{"D",6,3,0,-25},{"E",8,0,0,0} };
        std::vector<Barra> barras = { hacer_barra(nodos,"A","B"),hacer_barra(nodos,"B","C"),hacer_barra(nodos,"A","C"),hacer_barra(nodos,"B","D"),hacer_barra(nodos,"C","D"),hacer_barra(nodos,"D","E"),hacer_barra(nodos,"C","E") };
        std::vector<Apoyo> apoyos = { {idx_nodo(nodos,"A"),TipoApoyo::PASADOR},{idx_nodo(nodos,"E"),TipoApoyo::RODILLO_Y} };
        auto res = resolver_armadura(nodos, barras, apoyos);
        imprimir_resultado("AC.2 | Warren PB=15 PD=25 kN", res, nodos);
        exportar_csv("AC2_resultados.csv", "AC.2", res, nodos);
    }
    // AC3 Px=8 Py=20
    {
        std::vector<Nodo> nodos = { {"A",0,0,0,0},{"B",3,4,8,-20},{"C",6,0,0,0},{"D",3,0,0,0} };
        std::vector<Barra> barras = { hacer_barra(nodos,"A","B"),hacer_barra(nodos,"B","D"),hacer_barra(nodos,"B","C"),hacer_barra(nodos,"A","D"),hacer_barra(nodos,"D","C") };
        std::vector<Apoyo> apoyos = { {idx_nodo(nodos,"A"),TipoApoyo::PASADOR},{idx_nodo(nodos,"C"),TipoApoyo::RODILLO_Y} };
        auto res = resolver_armadura(nodos, barras, apoyos);
        imprimir_resultado("AC.3 | Oblicua Px=8 Py=20 kN", res, nodos);
        exportar_csv("AC3_resultados.csv", "AC.3", res, nodos);
    }
    // AC4 P=10
    {
        std::vector<Nodo> nodos = { {"A",0,0,0,0},{"C",3,0,0,0},{"E",6,0,0,0},{"G",9,0,0,0},{"I",12,0,0,0},{"B",3,3,0,-10},{"D",6,3,0,-10},{"F",9,3,0,-10} };
        std::vector<Barra> barras = { hacer_barra(nodos,"A","C"),hacer_barra(nodos,"C","E"),hacer_barra(nodos,"E","G"),hacer_barra(nodos,"G","I"),hacer_barra(nodos,"B","D"),hacer_barra(nodos,"D","F"),hacer_barra(nodos,"C","B"),hacer_barra(nodos,"E","D"),hacer_barra(nodos,"G","F"),hacer_barra(nodos,"A","B"),hacer_barra(nodos,"F","I"),hacer_barra(nodos,"B","E"),hacer_barra(nodos,"D","G") };
        std::vector<Apoyo> apoyos = { {idx_nodo(nodos,"A"),TipoApoyo::PASADOR},{idx_nodo(nodos,"I"),TipoApoyo::RODILLO_Y} };
        auto res = resolver_armadura(nodos, barras, apoyos);
        imprimir_resultado("AC.4 | Pratt P=10 kN", res, nodos);
        exportar_csv("AC4_resultados.csv", "AC.4 P=10kN", res, nodos);
    }
}

// ─────────────────────────────────────────────────────────────
//  MENU PRINCIPAL
// ─────────────────────────────────────────────────────────────

int main() {
    while (true) {
        std::cout << "\n============================================================\n";
        std::cout << "  ANALISIS DE ARMADURAS PLANAS - METODO DE NODOS\n";
        std::cout << "  Practica ICC | CETYS Universidad\n";
        std::cout << "============================================================\n";
        std::cout << "\n  1) Ingresar armadura manualmente\n";
        std::cout << "  2) Caso AC.1 - Triangular basica\n";
        std::cout << "  3) Caso AC.2 - Warren de dos claros\n";
        std::cout << "  4) Caso AC.3 - Carga oblicua resultante\n";
        std::cout << "  5) Caso AC.4 - Pratt + estudio parametrico\n";
        std::cout << "  6) Ejecutar todos los casos (valores predeterminados)\n";
        std::cout << "  0) Salir\n";
        std::cout << "\n  Opcion: ";

        int op;
        if (!(std::cin >> op)) {
            std::cin.clear();
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            continue;
        }

        try {
            switch (op) {
            case 0: std::cout << "\n  Saliendo...\n"; return 0;
            case 1: ingresar_armadura_manual(); break;
            case 2: caso_AC1(); break;
            case 3: caso_AC2(); break;
            case 4: caso_AC3(); break;
            case 5: caso_AC4(); break;
            case 6: ejecutar_todos(); break;
            default: std::cout << "  Opcion no valida.\n";
            }
        }
        catch (const std::exception& e) {
            std::cout << "\n  [ERROR] " << e.what() << "\n";
        }

        std::cout << "\n  Presione Enter para continuar...";
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        std::cin.get();
    }
}