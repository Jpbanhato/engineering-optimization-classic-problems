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

// -----------------------------------------------------------

// DEFS
#define POP_SIZE 50
#define NUM_OF_INDEPENDENT_RUNS 35
#define NUM_OF_VARIABLES 3
#define NUM_OF_RESTRICOES 4

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

    if (best_sol != nullptr) {
        delete best_sol;
    }

    if (f_pais != nullptr) {
        delete f_pais;
    }
    
    if (f_gerados != nullptr) {
        delete f_gerados;
    }

    if (prox_populacao != nullptr) {
        for (int i = 0; i < POP_SIZE; i++) 
            delete[] prox_populacao[i];
        delete[] prox_populacao;
    }

    //avaliacoes
    evals = 0;
    
    // depois instanciar os parametros referentes a nova run 
    best_sol = new double[NUM_OF_VARIABLES];    
    populacao = new double*[POP_SIZE];
    prox_populacao = new double*[POP_SIZE];
    fitness_individuos = new double[POP_SIZE];
    prox_fitness_individuos = new double[POP_SIZE];
    f_pais       = new double[POP_SIZE];
    f_gerados    = new double[POP_SIZE];
    viol_pais    = new double*[POP_SIZE];
    viol_gerados = new double*[POP_SIZE];
    for (int i = 0; i < POP_SIZE; i++) {
        viol_pais[i]    = new double[NUM_OF_RESTRICOES];
        viol_gerados[i] = new double[NUM_OF_RESTRICOES];
    }
    contador = -1;

    // parametros relacionados ao teste de reinicio qndo ha estagnacao
    ultimo_best = 1e10;
}


// DE:
#define CR 0.9
#define F 0.6

// ------------------------------------------------------------

// num avaliacoes funcao objetivo
int num_evals = 36000;

// variaveis do projeto
double N, D, d;
struct problema {
    std::string nome;
    int dimensao;
    double *lim_inf;
    double *lim_sup;
    double *ponteiro_F;
} problema;

// espirais ativos: N
double int_N_inf = 2.0;
double int_N_sup = 15.0;

// diametro volta: D
double int_D_inf = 0.25;
double int_D_sup = 1.3;

// diametro arame: d
double int_d_inf = 0.05;
double int_d_sup = 2;

// em vetor:
double *lim_inferior = new double[NUM_OF_VARIABLES]{int_N_inf, int_D_inf, int_d_inf};
double *lim_superior = new double[NUM_OF_VARIABLES]{int_N_sup, int_D_sup, int_d_sup};

// ------------------------------------------------------------

// volume
double V(double x1, double x2, double x3) {
    return (x1 + 2) * x2 * pow(x3, 2);
}

// restricoes mecanicas
double g1(double x1,double x2,double x3){ return 1.0 - (pow(x2,3.0)*x1)/(71785.0*pow(x3,4.0)); }

double g2(double x1,double x2,double x3){ return (4.0*x2*x2 - x3*x2)/(12566.0*(x2*pow(x3,3.0)-pow(x3,4.0))) + 1.0/(5108.0*x3*x3) - 1.0; }

double g3(double x1,double x2,double x3){ return 1.0 - (140.45*x3)/(x2*x2*x1); }

double g4(double x1,double x2,double x3){ return (x2 + x3)/1.5 - 1.0; }

// ------------------------------------------------------------

// gera pop inicial: [N, D, d]
void preenche_solucao_candidata(double *solucao_candidata) {
    for (int i = 0; i < NUM_OF_VARIABLES; i++) {
        solucao_candidata[i] = lim_inferior[i] + (lim_superior[i]-lim_inferior[i]) * ((double)rand()/RAND_MAX);
    }
}

void inicializa_populacao(double **populacao, double**prox_populacao) {
  // para cada solução candidata (lista com os nodes) da população (lista de listas)
  for (int i = 0; i < POP_SIZE; i++) {
      // preencher a lista interna (solucao candidata) com os nós de maneira aleatória
      double *solucao_candidata = new double[NUM_OF_VARIABLES];
      populacao[i] = solucao_candidata;
      preenche_solucao_candidata(solucao_candidata);
      // instancia sub vetores da prox populacao
      prox_populacao[i] = new double[NUM_OF_VARIABLES];
  }
}

double fitness_evaluation(double *solucao_candidata) {
    evals++; //conta avaliacao da funcao objetivo
    return V(solucao_candidata[0], solucao_candidata[1], solucao_candidata[2]);
}

// avalia funcao objetivo e as violacoes 
void avalia_mede_violacoes(double **pop, double *f_out, double **viol_out) {
    for (int i = 0; i < POP_SIZE; i++) {
        double *x = pop[i];
        f_out[i] = fitness_evaluation(x);
        viol_out[i][0] = std::max(0.0, g1(x[0],x[1],x[2]));
        viol_out[i][1] = std::max(0.0, g2(x[0],x[1],x[2]));
        viol_out[i][2] = std::max(0.0, g3(x[0],x[1],x[2]));
        viol_out[i][3] = std::max(0.0, g4(x[0],x[1],x[2]));
    }
}

// estatísticas dos pais
void calc_APM_stats(double *f_ref, double **viol_ref, double &media_f, double *media_v, double &soma_q) {
    media_f = 0.0;
    for (int j = 0; j < NUM_OF_RESTRICOES; j++)
        media_v[j] = 0.0;
    for (int i = 0; i < POP_SIZE; i++) {
        media_f += f_ref[i];
        for (int j = 0; j < NUM_OF_RESTRICOES; j++)
            media_v[j] += viol_ref[i][j];
    }
    media_f /= POP_SIZE;
    for (int j = 0; j < NUM_OF_RESTRICOES; j++)
        media_v[j] /= POP_SIZE;
    soma_q = 0.0;
    for (int j = 0; j < NUM_OF_RESTRICOES; j++)
        soma_q += media_v[j]*media_v[j];
}

// penaliza um indivíduo com as stats dadas
double penalizar(double f_i, double *viol_i, double media_f, double *media_v, double soma_q) {
    bool factivel = true;
    for (int j = 0; j < NUM_OF_RESTRICOES; j++)
        if (viol_i[j] > 0.0) {
            factivel = false; 
            break; 
        }
    if (factivel) 
        return f_i;
    double f_ = (f_i > media_f) ? f_i : media_f;
    double soma = 0.0;
    if (soma_q > 0.0)
        for (int j = 0; j < NUM_OF_RESTRICOES; j++)
            soma += (std::abs(media_f) * media_v[j] / soma_q) * viol_i[j];
    return f_ + soma;
}

int *sorteia_3_individuos(double **populacao, int i) {
    int *posicoes = new int[NUM_OF_VARIABLES];  
    do { posicoes[0] = rand() % POP_SIZE; } while (posicoes[0] == i);
    do { posicoes[1] = rand() % POP_SIZE; } while (posicoes[1] == posicoes[0] || posicoes[1] == i);
    do { posicoes[2] = rand() % POP_SIZE; } while (posicoes[2] == posicoes[0] || posicoes[2] == posicoes[1] || posicoes[2] == i);
    return posicoes;
}

// algoritmo do DE
void DE() {
    inicializa_populacao(populacao, prox_populacao);
    avalia_mede_violacoes(populacao, f_pais, viol_pais);

    double media_f, soma_q, *media_v = new double[NUM_OF_RESTRICOES];
    double melhor_run = 1e9;

    // melhor factivel da população inicial
    for (int i = 0; i < POP_SIZE; i++) {
        bool factivel = true;
        for (int j = 0; j < NUM_OF_RESTRICOES; j++)
            if (viol_pais[i][j] > 0.0)
                factivel = false;
        if (factivel && f_pais[i] < melhor_run) {
            melhor_run = f_pais[i];
            for (int j = 0; j < NUM_OF_VARIABLES; j++)
                best_sol[j] = populacao[i][j];
        }
    }

    int geracao = 0;
    while (evals < num_evals) {
        geracao++;
        // 1) gera os individuos novos
        for (int i = 0; i < POP_SIZE; i++) {
            int *idx = sorteia_3_individuos(populacao, i);
            int jrand = rand() % NUM_OF_VARIABLES;
            for (int j = 0; j < NUM_OF_VARIABLES; j++) {
                // crossover
                if (((double)rand()/RAND_MAX) < CR || j == jrand)
                    prox_populacao[i][j] = populacao[idx[2]][j] + F*(populacao[idx[0]][j] - populacao[idx[1]][j]);
                else
                    prox_populacao[i][j] = populacao[i][j];
                // ajusta os limites superior e inferior
                if (prox_populacao[i][j] < lim_inferior[j])
                    prox_populacao[i][j] = lim_inferior[j];
                if (prox_populacao[i][j] > lim_superior[j])
                    prox_populacao[i][j] = lim_superior[j];
            }
            delete[] idx;
        }

        // 2) mede para os individuos gerados
        avalia_mede_violacoes(prox_populacao, f_gerados, viol_gerados);

        // 3) mede estatisticas do APM sobre os pais
        calc_APM_stats(f_pais, viol_pais, media_f, media_v, soma_q);

        // 4) seleção com penalização consistente
        for (int i = 0; i < POP_SIZE; i++) {
            double F_pai = penalizar(f_pais[i], viol_pais[i], media_f, media_v, soma_q);
            double F_gerado = penalizar(f_gerados[i], viol_gerados[i], media_f, media_v, soma_q);
            if (F_gerado <= F_pai) {
                for (int j = 0; j < NUM_OF_VARIABLES; j++)
                    populacao[i][j] = prox_populacao[i][j];
                f_pais[i] = f_gerados[i];
                for (int j = 0; j < NUM_OF_RESTRICOES; j++)
                    viol_pais[i][j] = viol_gerados[i][j];

                // arquiva melhor factivel
                bool factivel = true;
                for (int j = 0; j < NUM_OF_RESTRICOES; j++)
                    if (viol_pais[i][j] > 0.0)
                        factivel = false;
                if (factivel && f_pais[i] < melhor_run) {
                    melhor_run = f_pais[i];
                    for (int j = 0; j < NUM_OF_VARIABLES; j++)
                        best_sol[j] = populacao[i][j];
                }
            }
        }
        // if (geracao % 100 == 0) {
        //     int nfact = 0;
        //     for (int i = 0; i < POP_SIZE; i++) {
        //         bool f = true;
        //         for (int j = 0; j < NUM_OF_RESTRICOES; j++) if (viol_pais[i][j] > 0.0) f = false;
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

//main
int main() {
    double resultados[NUM_OF_INDEPENDENT_RUNS];
    double best_global = 1e30, best_vars[NUM_OF_VARIABLES];
    for (int run = 0; run < NUM_OF_INDEPENDENT_RUNS; run++) {
        srand(run);
        inicializa_DE();
        DE();
        resultados[run] = ultimo_best;
        if (ultimo_best < best_global) {
        best_global = ultimo_best;
        for (int j = 0; j < NUM_OF_VARIABLES; j++) 
            best_vars[j] = best_sol[j];
    }
        printf("run %2d: %.6f\n", run+1, ultimo_best);
        printf("vars:  N=%.5f  D=%.5f  d=%.5f  (V=%.6f)\n", best_vars[0], best_vars[1], best_vars[2], best_global);
    }
    std::sort(resultados, resultados+NUM_OF_INDEPENDENT_RUNS);
    double soma = 0;
    for (double r : resultados)
        soma += r;
    double media = soma/NUM_OF_INDEPENDENT_RUNS, dp = 0;
    for (double r : resultados)
        dp += (r-media)*(r-media);
    dp = sqrt(dp/NUM_OF_INDEPENDENT_RUNS);
    printf("melhor=%.6f mediana=%.6f media=%.6f dp=%.4e pior=%.6f\n", resultados[0], resultados[17], media, dp, resultados[34]);
    printf("vars:  N=%.5f  D=%.5f  d=%.5f  (V=%.6f)\n", best_vars[0], best_vars[1], best_vars[2], best_global);
    return 0;
}








