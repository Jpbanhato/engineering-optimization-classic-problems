#include<cmath>
#include<iostream>
#include<algorithm>
#include<stdio.h>
#include<stdlib.h>
#include<string>
#include<cstring>
#include<math.h>
#include<fstream>
#include<limits.h>
#include<string>
#include<stdexcept>
#include<unordered_map>
#include<vector>

// -----------------------------------------------------------

// DEFS
#define POP_SIZE 50
#define NUM_OF_INDEPENDENT_RUNS 35

// DE:
#define CR 0.9
#define F 0.6


// PROBLEMA
struct problema {
    std::string nome;
    int num_variaveis;
    const char *variaveis[10];
    int num_restricoes;
    std::vector<double(*)(double*)> restricoes;
    double *lim_inf;
    double *lim_sup;
    double (*funcao_objetivo)(double*);
    bool *continuas;
    double *passo;
    int num_evals;
} problema;

// -----------------------------------------------------------


// gerais
static double **populacao = nullptr;
static double **prox_populacao = nullptr;
static int contador;
static double *fitness_individuos = nullptr;
static double *prox_fitness_individuos = nullptr;
static double ultimo_best;
double *best_sol = nullptr;
static int evals = 0;

static double  *f_pais        = nullptr;
static double **viol_pais           = nullptr;
static double  *f_gerados  = nullptr;
static double **viol_gerados     = nullptr;

static double  *xproj = nullptr;

void inicializa_DE() {
    // primeiro limpar os vetores dinamicos
    if (populacao != nullptr) {
        for (int i = 0 ; i < POP_SIZE; i++)
            delete[] populacao[i];
        delete[] populacao;
    }
    if (viol_pais != nullptr) {
        for (int i = 0 ; i < POP_SIZE; i++)
            delete[] viol_pais[i];
        delete[] viol_pais;
    }
    if (viol_gerados != nullptr) {
        for (int i = 0 ; i < POP_SIZE; i++)
            delete[] viol_gerados[i];
        delete[] viol_gerados;
    }

    if (fitness_individuos != nullptr)
        delete[] fitness_individuos;

    if (prox_fitness_individuos != nullptr) 
        delete[] prox_fitness_individuos;

    if (best_sol != nullptr) { delete[] best_sol; }
    if (f_pais != nullptr)   { delete[] f_pais; }
    if (f_gerados != nullptr){ delete[] f_gerados; }
    
    if (xproj != nullptr)   { delete[] xproj; }
    
    if (prox_populacao != nullptr) {
        for (int i = 0; i < POP_SIZE; i++) 
            delete[] prox_populacao[i];
        delete[] prox_populacao;
    }

    //avaliacoes
    evals = 0;
    
    // depois instanciar os parametros referentes a nova run 
    best_sol = new double[problema.num_variaveis];    
    populacao = new double*[POP_SIZE];
    prox_populacao = new double*[POP_SIZE];
    fitness_individuos = new double[POP_SIZE];
    prox_fitness_individuos = new double[POP_SIZE];
    f_pais       = new double[POP_SIZE];
    
    xproj       = new double[problema.num_variaveis];

    f_gerados    = new double[POP_SIZE];
    viol_pais    = new double*[POP_SIZE];
    viol_gerados = new double*[POP_SIZE];
    for (int i = 0; i < POP_SIZE; i++) {
        viol_pais[i]    = new double[problema.num_restricoes];
        viol_gerados[i] = new double[problema.num_restricoes];
    }
    contador = -1;

    // parametros relacionados ao teste de reinicio qndo ha estagnacao
    ultimo_best = 1e10;
}

// ------------------------------------------------------------

// gera pop inicial: [x1, x2, ..., xn]
void preenche_solucao_candidata(double *solucao_candidata) {
    for (int i = 0; i < problema.num_variaveis; i++) {
        solucao_candidata[i] = problema.lim_inf[i] + (problema.lim_sup[i]-problema.lim_inf[i]) * ((double)rand()/RAND_MAX);
    }
}

void inicializa_populacao(double **populacao, double**prox_populacao) {
  for (int i = 0; i < POP_SIZE; i++) {
      double *solucao_candidata = new double[problema.num_variaveis];
      populacao[i] = solucao_candidata;
      preenche_solucao_candidata(solucao_candidata);
      prox_populacao[i] = new double[problema.num_variaveis];
  }
}

double fitness_evaluation(double *solucao_candidata) {
    evals++; //conta avaliacao da funcao objetivo
    return problema.funcao_objetivo(solucao_candidata);
}

// para tratar de problema misto (variavel se torna um multiplo do passo, ou se multiplo for 0, se torna inteira)
void projeta(const double* indiv, double* indiv_projetado) {
    for (int j = 0; j < problema.num_variaveis; j++) {
        double v = indiv[j];
        if (problema.continuas[j] == false) { //ou seja, se a variavel é inteira / contem passo:
            if (problema.passo[j] > 0.0) //contem passo
                v = std::round(v / problema.passo[j]) * problema.passo[j];
            else //variavel inteira
                v = std::round(v);
        }
        indiv_projetado[j] = v;
    }
}

// avalia funcao objetivo e as violacoes 
void avalia_mede_violacoes(double **pop, double *f_out, double **viol_out) {
    for (int i = 0; i < POP_SIZE; i++) {
        double *x = pop[i];
        projeta(x, xproj);
        f_out[i] = fitness_evaluation(xproj); //funcao objetivo
        for (int v = 0; v < problema.num_restricoes; v++) {
            viol_out[i][v] = std::max(0.0, problema.restricoes[v](xproj)); //violacoes
        }
    }
}

// estatísticas dos pais
void calc_APM_stats(double *f_ref, double **viol_ref, double &media_f, double *media_v, double &soma_q) {
    media_f = 0.0;
    for (int j = 0; j < problema.num_restricoes; j++)
        media_v[j] = 0.0;
    for (int i = 0; i < POP_SIZE; i++) {
        media_f += f_ref[i];
        for (int j = 0; j < problema.num_restricoes; j++)
            media_v[j] += viol_ref[i][j];
    }
    media_f /= POP_SIZE;
    for (int j = 0; j < problema.num_restricoes; j++)
        media_v[j] /= POP_SIZE;
    soma_q = 0.0;
    for (int j = 0; j < problema.num_restricoes; j++)
        soma_q += media_v[j]*media_v[j];
}

// penaliza um indivíduo com as stats dadas
double penalizar(double f_i, double *viol_i, double media_f, double *media_v, double soma_q) {
    bool factivel = true;
    for (int j = 0; j < problema.num_restricoes; j++)
        if (viol_i[j] > 0.0) {
            factivel = false; 
            break; 
        }
    if (factivel) //se é factivel
        return f_i;
    // se nao é factivel
    double f_ = (f_i > media_f) ? f_i : media_f;
    double soma = 0.0;
    if (soma_q > 0.0)
        for (int j = 0; j < problema.num_restricoes; j++)
            soma += (std::abs(media_f) * media_v[j] / soma_q) * viol_i[j];
    return f_ + soma;
}

int *sorteia_3_individuos(double **populacao, int i) {
    int *posicoes = new int[3];  
    do { posicoes[0] = rand() % POP_SIZE; } while (posicoes[0] == i);
    do { posicoes[1] = rand() % POP_SIZE; } while (posicoes[1] == posicoes[0] || posicoes[1] == i);
    do { posicoes[2] = rand() % POP_SIZE; } while (posicoes[2] == posicoes[0] || posicoes[2] == posicoes[1] || posicoes[2] == i);
    return posicoes;
}

// algoritmo do DE
void DE() {
    inicializa_populacao(populacao, prox_populacao);
    avalia_mede_violacoes(populacao, f_pais, viol_pais);

    double media_f, soma_q, *media_v = new double[problema.num_restricoes];
    double melhor_run = 1e9;

    // melhor factivel da população inicial
    for (int i = 0; i < POP_SIZE; i++) {
        bool factivel = true;
        for (int j = 0; j < problema.num_restricoes; j++)
            if (viol_pais[i][j] > 0.0)
                factivel = false;
        if (factivel && f_pais[i] < melhor_run) {
            melhor_run = f_pais[i];
            for (int j = 0; j < problema.num_variaveis; j++)
                best_sol[j] = populacao[i][j];
        }
    }

    int geracao = 0;
    while (evals < problema.num_evals) {
        geracao++;
        // 1) gera os individuos novos
        for (int i = 0; i < POP_SIZE; i++) {
            int *idx = sorteia_3_individuos(populacao, i);
            int jrand = rand() % problema.num_variaveis;
            for (int j = 0; j < problema.num_variaveis; j++) {
                // crossover
                if (((double)rand()/RAND_MAX) < CR || j == jrand)
                    prox_populacao[i][j] = populacao[idx[2]][j] + F*(populacao[idx[0]][j] - populacao[idx[1]][j]);
                else
                    prox_populacao[i][j] = populacao[i][j];
                // ajusta os limites superior e inferior
                if (prox_populacao[i][j] < problema.lim_inf[j])
                    prox_populacao[i][j] = problema.lim_inf[j];
                if (prox_populacao[i][j] > problema.lim_sup[j])
                    prox_populacao[i][j] = problema.lim_sup[j];
            }
            delete[] idx;
        }

        // 2) mede para os individuos gerados
        avalia_mede_violacoes(prox_populacao, f_gerados, viol_gerados);

        // 3) mede estatisticas do APM sobre os pais: media_fitness_geracao_dos_pais, media_violacao_j_geracao_dos_pais, media_violacoes^2
        calc_APM_stats(f_pais, viol_pais, media_f, media_v, soma_q);

        // 4) seleção com penalização APM
        for (int i = 0; i < POP_SIZE; i++) {
            double F_pai = penalizar(f_pais[i], viol_pais[i], media_f, media_v, soma_q);
            double F_gerado = penalizar(f_gerados[i], viol_gerados[i], media_f, media_v, soma_q);
            if (F_gerado <= F_pai) {
                for (int j = 0; j < problema.num_variaveis; j++)
                    populacao[i][j] = prox_populacao[i][j];
                f_pais[i] = f_gerados[i];
                for (int j = 0; j < problema.num_restricoes; j++)
                    viol_pais[i][j] = viol_gerados[i][j];

                // arquiva melhor factivel
                bool factivel = true;
                for (int j = 0; j < problema.num_restricoes; j++)
                    if (viol_pais[i][j] > 0.0)
                        factivel = false;
                if (factivel && f_pais[i] < melhor_run) {
                    melhor_run = f_pais[i];
                    for (int j = 0; j < problema.num_variaveis; j++)
                        best_sol[j] = populacao[i][j];
                }
            }
        }
        // if (geracao % 100 == 0) {
        //     int nfact = 0;
        //     for (int i = 0; i < POP_SIZE; i++) {
        //         bool f = true;
        //         for (int j = 0; j < problema.num_restricoes; j++) if (viol_pais[i][j] > 0.0) f = false;
        //         if (f) nfact++;
        //     }
        //     fprintf(stderr, "  [gen %4d | evals %6d] melhor=%.6f | factiveis=%d/%d\n",
        //             geracao, evals, melhor_run, nfact, POP_SIZE);
        // }
    }
    delete[] media_v;
    ultimo_best = melhor_run;
}

// ------------------------------------------------------------
// ----------------PROBLEMA 1-------------------
// volume
double V_p1(double *x) { return (x[0] + 2) * x[1] * pow(x[2], 2); }

// restricoes mecanicas
double g1_p1(double *x){ return 1.0 - (pow(x[1],3.0)*x[0])/(71785.0*pow(x[2],4.0)); }

double g2_p1(double *x){ return (4.0*x[1]*x[1] - x[2]*x[1])/(12566.0*(x[1]*pow(x[2],3.0)-pow(x[2],4.0))) + 1.0/(5108.0*x[2]*x[2]) - 1.0; }

double g3_p1(double *x){ return 1.0 - (140.45*x[2])/(x[1]*x[1]*x[0]); }

double g4_p1(double *x){ return (x[1] + x[2])/1.5 - 1.0; }

// ----------------PROBLEMA 2-------------------
double W_p2 (double* x){
    return 0.7854*x[0]*x[1]*x[1]*(3.3333*x[2]*x[2] + 14.9334*x[2] - 43.0934) - 1.508*x[0]*(x[5]*x[5] + x[6]*x[6]) + 7.4777*(std::pow(x[5],3.0) + std::pow(x[6],3.0)) + 0.7854*(x[3]*x[5]*x[5] + x[4]*x[6]*x[6]);
}
double g1_p2 ( double* x){ return 27.0/(x[0]*x[1]*x[1]*x[2]) - 1.0; }
double g2_p2 ( double* x){ return 397.5/(x[0]*x[1]*x[1]*x[2]*x[2]) - 1.0; }
double g3_p2 ( double* x){ return 1.93*std::pow(x[3],3.0)/(x[1]*x[2]*std::pow(x[5],4.0)) - 1.0; }
double g4_p2 ( double* x){ return 1.93*std::pow(x[4],3.0)/(x[1]*x[2]*std::pow(x[6],4.0)) - 1.0; }
double g5_p2 ( double* x){ double t=745.0*x[3]/(x[1]*x[2]); return std::sqrt(t*t+16.9e6)/(110.0*std::pow(x[5],3.0)) - 1.0; }
double g6_p2 ( double* x){ double t=745.0*x[4]/(x[1]*x[2]); return std::sqrt(t*t+157.5e6)/(85.0*std::pow(x[6],3.0)) - 1.0; }
double g7_p2 ( double* x){ return x[1]*x[2]/40.0 - 1.0; }
double g8_p2 ( double* x){ return 5.0*x[1]/x[0] - 1.0; }          // x1/x2 >= 5
double g9_p2 ( double* x){ return x[0]/(12.0*x[1]) - 1.0; }       // x1/x2 <= 12
double g10_p2( double* x){ return (1.5*x[5]+1.9)/x[3] - 1.0; }
double g11_p2( double* x){ return (1.1*x[6]+1.9)/x[4] - 1.0; }

// ----------------PROBLEMA 3-------------------
// f obj
double C_p3 (double *x) { return 1.10471 * x[0]*x[0] * x[1] + 0.04811 * x[2] * x[3] * (14.0 + x[1]); }
//aux
double alpha_p3 (double *x) { return sqrt( 0.25 * (x[1]*x[1] + (x[0] + x[2])*(x[0] + x[2])) ); }
double tau_linha_p3 (double *x) { return (6000.0) / (sqrt(2) * x[0] * x[1]) ; }
double tau_2linha_p3 (double *x) { return (6000.0 * (14.0 + 0.5*x[1]) * alpha_p3(x)) / (2.0 * 0.707 * x[0] * x[1] * ( (x[1]*x[1] / 12.0) + 0.25 * pow(x[0]+x[2], 2) )) ; }
double tau_p3 (double *x) { return sqrt( pow(tau_linha_p3(x), 2) + pow(tau_2linha_p3(x), 2) + (x[1]*tau_linha_p3(x)*tau_2linha_p3(x))/(alpha_p3(x)) ); }
double ro_p3(double *x) { return 64746.022 * (1 - 0.0282346*x[2]) * x[2]*pow(x[3], 3) ; }
double sigma_p3(double *x) { return 504000.0 / (x[2]*x[2]*x[3]) ; }
double delta_p3(double *x) { return (2.1952) / (pow(x[2], 3) * x[3]) ; }
//rest
double g1_p3 (double *x) { return - (13600 - tau_p3(x)); }
double g2_p3 (double *x) { return - (30000 - sigma_p3(x)); }
double g3_p3 (double *x) { return - (x[3] - x[0]); }
double g4_p3 (double *x) { return - (ro_p3(x) - 6000); }
double g5_p3 (double *x) { return - (0.25 - delta_p3(x)); }

// ----------------PROBLEMA 3-------------------
double W_p4(double *x) { return 0.6224 * x[0] * x[2] * x[3] + 1.7781 * x[1] * x[2] * x[2] + 3.1661 * x[0] * x[0] * x[3] + 19.84  * x[0] * x[0] * x[2]; ;}
double g1_p4 (double *x) { return -(x[0] - 0.0193*x[2]) ; }
double g2_p4 (double *x) { return -(x[1] - 0.00954*x[2]) ; }
double g3_p4 (double *x) { return -((M_PI * x[2]*x[2]*x[3] + (4.0/3.0) * M_PI * x[2]*x[2]*x[2])/1296000.0 - 1.0) ; }
double g4_p4 (double *x) { return -(-x[3] + 240) ; }

// ---------------------------------------------
// ------------------------------------------------------------

void definir_problema(int num_problema) {
    if (num_problema == 1) {
        problema.nome = "mola sob tracao / compressao";
        
        problema.num_variaveis = 3;
        
        problema.variaveis[0] = "N";
        problema.variaveis[1] = "D";
        problema.variaveis[2] = "d";
        
        problema.num_restricoes = 4;
        problema.restricoes = { g1_p1, g2_p1, g3_p1, g4_p1 };
        problema.funcao_objetivo = V_p1;
        
        problema.lim_inf = new double[problema.num_variaveis];
        problema.lim_sup = new double[problema.num_variaveis];
        problema.lim_inf[0] = 2.0; problema.lim_sup[0] = 15.0;
        problema.lim_inf[1] = 0.25; problema.lim_sup[1] = 1.3;
        problema.lim_inf[2] = 0.05; problema.lim_sup[2] = 2.0;

        problema.continuas = new bool[problema.num_variaveis] {true, true, true};        
        problema.passo = new double[problema.num_variaveis] {0, 0, 0};
        problema.num_evals = 36000;
        
    } else if (num_problema == 2) {
        problema.nome = "redutor de velocidade";
        
        problema.num_variaveis = 7;
        
        problema.variaveis[0] = "b";
        problema.variaveis[1] = "m";
        problema.variaveis[2] = "n"; //inteira
        problema.variaveis[3] = "l1";
        problema.variaveis[4] = "l2";
        problema.variaveis[5] = "d1";
        problema.variaveis[6] = "d2";
        
        problema.num_restricoes = 11;
        problema.restricoes = { g1_p2, g2_p2, g3_p2, g4_p2, g5_p2, g6_p2, g7_p2, g8_p2, g9_p2, g10_p2, g11_p2 };
        problema.funcao_objetivo = W_p2;
        
        problema.lim_inf = new double[problema.num_variaveis];
        problema.lim_sup = new double[problema.num_variaveis];
        problema.lim_inf[0] = 2.6; problema.lim_sup[0] = 3.6;
        problema.lim_inf[1] = 0.7; problema.lim_sup[1] = 0.8;
        problema.lim_inf[2] = 17; problema.lim_sup[2] = 28; //inteira
        problema.lim_inf[3] = 7.3; problema.lim_sup[3] = 8.3;
        problema.lim_inf[4] = 7.8; problema.lim_sup[4] = 8.3;
        problema.lim_inf[5] = 2.9; problema.lim_sup[5] = 3.9;
        problema.lim_inf[6] = 5.0; problema.lim_sup[6] = 5.9;

        problema.continuas = new bool[problema.num_variaveis] {true, true, false, true, true, true, true};
        problema.passo = new double[problema.num_variaveis] {0, 0, 0, 0, 0, 0, 0};
        problema.num_evals = 36000;

    } else if (num_problema == 3) {
        problema.nome = "viga soldada";
        
        problema.num_variaveis = 4;
        
        problema.variaveis[0] = "h";
        problema.variaveis[1] = "l";
        problema.variaveis[2] = "t";
        problema.variaveis[3] = "b";
        
        problema.num_restricoes = 5;
        problema.restricoes = { g1_p3, g2_p3, g3_p3, g4_p3, g5_p3 };
        problema.funcao_objetivo = C_p3;
        
        problema.lim_inf = new double[problema.num_variaveis];
        problema.lim_sup = new double[problema.num_variaveis];
        problema.lim_inf[0] = 0.125; problema.lim_sup[0] = 10;
        problema.lim_inf[1] = 0.1; problema.lim_sup[1] = 10;
        problema.lim_inf[2] = 0.1; problema.lim_sup[2] = 10;
        problema.lim_inf[3] = 0.1; problema.lim_sup[3] = 10;

        problema.continuas = new bool[problema.num_variaveis] {true, true, true, true};
        problema.passo = new double[problema.num_variaveis] {0, 0, 0, 0};
        problema.num_evals = 320000;

    } else if (num_problema == 4) {
        problema.nome = "vaso de pressao";
        
        problema.num_variaveis = 4;
        
        problema.variaveis[0] = "Ts"; //inteiro
        problema.variaveis[1] = "Th"; //inteiro
        problema.variaveis[2] = "R";
        problema.variaveis[3] = "L";
        
        problema.num_restricoes = 4;
        problema.restricoes = { g1_p4, g2_p4, g3_p4, g4_p4 };
        problema.funcao_objetivo = W_p4;
        
        problema.lim_inf = new double[problema.num_variaveis];
        problema.lim_sup = new double[problema.num_variaveis];
        problema.lim_inf[0] = 0.00625; problema.lim_sup[0] = 5;
        problema.lim_inf[1] = 0.00625; problema.lim_sup[1] = 5;
        problema.lim_inf[2] = 10; problema.lim_sup[2] = 200;
        problema.lim_inf[3] = 10; problema.lim_sup[3] = 200;

        problema.continuas = new bool[problema.num_variaveis] {false, false, true, true};
        problema.passo = new double[problema.num_variaveis] {0.0625, 0.0625, 0, 0};
        problema.num_evals = 80000;
    } else {
        throw std::runtime_error("Indice do problema invalido.");
    }
}

//main
int main() {
    // problema
    int num_problema = 4;
    definir_problema(num_problema);
    printf("===========================================================================\n");
    printf("SOLUCIONANDO PROBLEMA: %s \n", problema.nome.c_str());
    printf("===========================================================================\n");
    // solucao
    double resultados[NUM_OF_INDEPENDENT_RUNS];
    double best_global = 1e30, best_vars[problema.num_variaveis];
    for (int run = 0; run < NUM_OF_INDEPENDENT_RUNS; run++) {
        srand(run);
        inicializa_DE();
        DE();
        resultados[run] = ultimo_best;
        if (ultimo_best < best_global) {
            best_global = ultimo_best;
            for (int j = 0; j < problema.num_variaveis; j++) 
                best_vars[j] = best_sol[j];
        }
        printf("run %2d: %.6f\n", run+1, ultimo_best);
        printf("vars: ");
        if (ultimo_best < best_global) {
            best_global = ultimo_best;
            projeta(best_sol, best_vars);
        }
        for (int i = 0; i < problema.num_variaveis; i++)
            printf(" %s=%.5f ", problema.variaveis[i], best_vars[i]);
        printf(" (best_global=%.6f)", best_global);
        printf("\n");
    }
    std::sort(resultados, resultados+NUM_OF_INDEPENDENT_RUNS);
    double soma = 0;
    for (double r : resultados)
        soma += r;
    double media = soma/NUM_OF_INDEPENDENT_RUNS, dp = 0;
    for (double r : resultados)
        dp += (r-media)*(r-media);
    dp = sqrt(dp/NUM_OF_INDEPENDENT_RUNS);
    printf("-----> melhor=%.6f mediana=%.6f media=%.6f dp=%.4e pior=%.6f\n", resultados[0], (NUM_OF_INDEPENDENT_RUNS % 2 != 0) ? resultados[NUM_OF_INDEPENDENT_RUNS / 2] : (resultados[NUM_OF_INDEPENDENT_RUNS / 2 - 1] + resultados[NUM_OF_INDEPENDENT_RUNS / 2]) / 2.0, media, dp, resultados[NUM_OF_INDEPENDENT_RUNS - 1]);
    
    printf("vars: ");
    for (int i = 0; i < problema.num_variaveis; i++)
        printf(" %s=%.5f ", problema.variaveis[i], best_vars[i]);
    printf(" (best_global=%.6f)", best_global);
    printf("\n");
    return 0;
}


