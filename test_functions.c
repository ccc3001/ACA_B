#include "new_aca.h"
#include "aca.h"
#include "aca_b.h"
#include "fullmatrix.h"
#include "rkmatrix.h"
#include "basic.h"
#include "kernel_functions.h"
#include "svd.h"
#include "interpolation.h"
#include "svd.h"
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>
#include <cjson/cJSON.h>
#include <string.h>
//   input matrix
//   output matrix
//   timestamp
//   residuals

#define JSON_STRING(obj, key) \
    (cJSON_GetObjectItem(obj, key)->valuestring)

int matrix_aca_test(int d, int L, bool runtime_benchmark, const char *test_function_name, int iter, int repeats, double eps,
                    bool residuals_benchmark, double eps_residual, int starts_res, int dim, int rows, int cols, const double *nodes_x, const double *nodes_y,
                    double (*test_function)(int dim, const double *x, const double *y),bool hbaca_residual_test,double eps_hbaca,int h_baca_repeats)
{
    time_t now = time(NULL);
    struct tm *t = localtime(&now);

    char filename[100];

    strftime(filename, sizeof(filename), "results/%Y-%m-%d_%H-%M-%S.json", t);

    cJSON *root = cJSON_CreateObject();
    if (d == 0 && L == 0)
    {
        cJSON_AddStringToObject(root, "method", "ACA");
    }
    else if (d > 0 && L == 0)
    {
        cJSON_AddStringToObject(root, "method", "BACA");
    }
    else if (d > 0 && L > 0)
    {
        cJSON_AddStringToObject(root, "method", "HBACA");
    }
    else
    {
        printf("Error when choosing method");
        cJSON_Delete(root);
        return 0;
    }

    cJSON_AddStringToObject(root, "test_function", test_function_name);
    cJSON_AddStringToObject(root, "timestamp", filename);

    // cJSON_AddNumberToObject(root, "testfunc", testfunc);
    cJSON_AddNumberToObject(root, "d", d);
    cJSON_AddNumberToObject(root, "L", L);
    if (runtime_benchmark)
    {
        cJSON *bench = cJSON_CreateObject();

        cJSON_AddBoolToObject(bench, "enabled", true);
        cJSON_AddNumberToObject(bench, "epsilon", eps);

        cJSON_AddNumberToObject(bench, "iterations", iter);
        cJSON_AddNumberToObject(bench, "repeats", repeats);

        cJSON *cases = cJSON_CreateArray();
        cJSON_AddItemToObject(bench, "cases", cases);

        printf("\n-------------------------------------\n");
        printf("Running %s%s", JSON_STRING(root, "method"), "\n");
        for (int i = 0; i < iter; i++)
        {
            cJSON *case_json = cJSON_CreateObject();

            cJSON *runs = cJSON_CreateArray();
            int n = 1000 + (1 + i) * 100;
            int dim_ = 1;
            int nd_test = n;
            double *nodes_x_test = allocate_doubles(dim_ * nd_test);
            double *nodes_y_test = allocate_doubles(dim_ * nd_test);

            double h_ = 1.0 / (n + 1);

            for (int k = 0; k < nd_test; k++)
            {
                nodes_x_test[k] = (k + 1) * h_;

                // nodes_x_test[nd_test + k] = (k + 1) * h_;

                nodes_y_test[k] = 0.2 + (k + 1) * h_;

                // nodes_y_test[nd_test + k] = 0.5 + (k + 1) * h_;
            }
            printf("test");
            cJSON_AddNumberToObject(case_json, "n", n);
            double total_time = 0.0;
            double total_rank = 0.0;

            double total_T_over_rank = 0.0;
            double total_T_over_nk = 0.0;
            double total_T_over_nk2 = 0.0;
            printf("\n=====================================\n");
            printf("Matrix size: %d x %d\n", n, n);
            printf("=====================================\n");

            for (int r = 0; r < repeats; r++)
            {
                printf("\nRepeat %d / %d\n", r + 1, repeats);
                cJSON *run = cJSON_CreateObject();

                prkmatrix RK_ACA = NULL;
                pfullmatrix U_h = NULL;
                pfullmatrix V_h = NULL;
                pfullmatrix S_h = NULL;

                double *residuals_u = NULL;
                double *residuals_v = NULL;
                double *rank_increase = NULL;

                double end, start;

                int rank = 0;

                if (strcmp(JSON_STRING(root, "method"), "ACA") == 0)
                {
                    start = omp_get_wtime();
                    RK_ACA = aca_rkmatrix_(eps, dim_, nd_test, nd_test, nodes_x_test, nodes_y_test, test_function, &residuals_u, &residuals_v);
                    cJSON_AddNumberToObject(case_json, "rows", n);
                    cJSON_AddNumberToObject(case_json, "cols", n);
                    cJSON_AddNumberToObject(case_json, "dim", dim_);
                    end = omp_get_wtime();
                    // TODO: this is not a good solution we need to calculate nodes_x and nodes_y not use the ones given by the function
                    cJSON *u_json = cJSON_CreateDoubleArray(residuals_u, RK_ACA->kt);
                    cJSON *v_json = cJSON_CreateDoubleArray(residuals_v, RK_ACA->kt);
                    cJSON_AddItemToObject(run, "residuals_u", u_json);
                    cJSON_AddItemToObject(run, "residuals_v", v_json);
                    rank = RK_ACA->kt;
                }
                else if (strcmp(JSON_STRING(root, "method"), "BACA") == 0)
                {
                    start = omp_get_wtime();
                    RK_ACA = b_aca_rkmatrix(eps, d, dim_, nd_test, nd_test, nodes_x_test, nodes_y_test, test_function, &residuals_u, &residuals_v, &rank_increase);
                    end = omp_get_wtime();
                    cJSON_AddNumberToObject(case_json, "rows", n);
                    cJSON_AddNumberToObject(case_json, "cols", n);
                    cJSON_AddNumberToObject(case_json, "dim", dim_);
                    cJSON *u_json = cJSON_CreateDoubleArray(residuals_u, RK_ACA->kt);
                    cJSON *v_json = cJSON_CreateDoubleArray(residuals_v, RK_ACA->kt);
                    cJSON *rank_increase_json = cJSON_CreateDoubleArray(rank_increase, RK_ACA->kt);
                    cJSON_AddItemToObject(run, "rank_increase", rank_increase_json);
                    cJSON_AddItemToObject(run, "residuals_u", u_json);
                    cJSON_AddItemToObject(run, "residuals_v", v_json);
                    rank = RK_ACA->kt;
                }
                else if (strcmp(JSON_STRING(root, "method"), "HBACA") == 0)
                {

                    int r_h;
                    ACAResidualNode *new_residual_node(void);
                    ACAResidualNode *root = new_residual_node();
                    root->row_start = 0;
                    root->col_start = 0;
                    start = omp_get_wtime();
                    h_b_aca_rkmatrix(eps, d, L, root, &U_h, &S_h, &V_h, dim_, nd_test, nd_test, nodes_x_test, nodes_y_test, test_function, &r_h);
                    end = omp_get_wtime();
                    cJSON_AddNumberToObject(case_json, "rows", n);
                    cJSON_AddNumberToObject(case_json, "cols", n);
                    cJSON_AddNumberToObject(case_json, "dim", dim_);
                    cJSON *tree_json = aca_residual_node_to_json(root);
                    cJSON_AddItemToObject(run, "h_b_aca_tree", tree_json);
                    del_residual_node(root);
                    rank = r_h;
                }

                double elapsed = end - start;
                printf("Elapsed: %.6f s | Rank: %d\n",
                       elapsed,
                       rank);
                double T_over_rank = 0.0;
                double T_over_nk = 0.0;
                double T_over_nk2 = 0.0;

                if (rank > 0)
                {
                    T_over_rank = elapsed / rank;

                    T_over_nk = elapsed / (n * rank);

                    T_over_nk2 = elapsed / (n * rank * rank);
                }

                total_time += elapsed;
                total_rank += rank;

                total_T_over_rank += T_over_rank;
                total_T_over_nk += T_over_nk;
                total_T_over_nk2 += T_over_nk2;
                cJSON_AddNumberToObject(run, "repeat", r);
                cJSON_AddNumberToObject(run, "elapsed", elapsed);
                cJSON_AddNumberToObject(run, "rank", rank);
                cJSON_AddNumberToObject(run, "T_over_rank", T_over_rank);
                cJSON_AddNumberToObject(run, "T_over_nk", T_over_nk);
                cJSON_AddNumberToObject(run, "T_over_nk2", T_over_nk2);
                cJSON_AddItemToArray(runs, run);
                if (strcmp(JSON_STRING(root, "method"), "HBACA") == 0)
                {
                    del_fullmatrix(U_h);
                    del_fullmatrix(V_h);
                    del_fullmatrix(S_h);
                }
                else
                {
                    del_rkmatrix(RK_ACA);
                }
            }

            double avg_time = total_time / repeats;

            double avg_rank = total_rank / repeats;

            double avg_T_over_rank = total_T_over_rank / repeats;

            double avg_T_over_nk = total_T_over_nk / repeats;

            double avg_T_over_nk2 = total_T_over_nk2 / repeats;

            cJSON *agg = cJSON_CreateObject();

            cJSON_AddNumberToObject(agg, "avg_time", avg_time);

            cJSON_AddNumberToObject(agg, "avg_rank", avg_rank);

            cJSON_AddNumberToObject(agg, "avg_T_over_rank", avg_T_over_rank);

            cJSON_AddNumberToObject(agg, "avg_T_over_nk", avg_T_over_nk);

            cJSON_AddNumberToObject(agg, "avg_T_over_nk2", avg_T_over_nk2);
            cJSON_AddItemToObject(case_json, "runs", runs);
            cJSON_AddItemToObject(case_json, "aggregated_runs", agg);
            cJSON_AddItemToArray(cases, case_json);
            free(nodes_x_test);
            free(nodes_y_test);
        }
        cJSON_AddItemToObject(root, "runtime_benchmark", bench);
    }

    if (residuals_benchmark)
    {
        printf("\n-------------------------------------\n");
        printf("Calculating residuals for %s%s", JSON_STRING(root, "method"), "\n");
        cJSON *residual_bench = cJSON_CreateObject();
        cJSON *starts_json = cJSON_CreateObject();
        cJSON_AddBoolToObject(residual_bench, "enabled", true);
        cJSON_AddNumberToObject(residual_bench, "epsilon", eps_residual);
        double best_aca_final_norm = HUGE_VAL;
        double *best_aca_frob_norms = NULL;

        int best_aca_rank = 0;

        int best_aca_start = -1;
        bool stop_after_this_run = false;
        for (int start = 0; start < starts_res; start++)
        {
            printf("\nACA random start %d / %d\n", start + 1, starts_res);
            cJSON *start_json = cJSON_CreateObject();
            int seed = time(NULL) + start;
            srand(seed);
            cJSON_AddNumberToObject(start_json, "seed", seed);
            double *residuals_u = NULL;
            double *residuals_v = NULL;
            double *rank_increase = NULL;
            prkmatrix RK_ACA = NULL;

            pfullmatrix U_h = NULL;
            pfullmatrix V_h = NULL;
            pfullmatrix S_h = NULL;

            int rank = 0;
            if (strcmp(JSON_STRING(root, "method"), "ACA") == 0)
            {
                RK_ACA = aca_rkmatrix_(eps_residual, dim, rows, cols, nodes_x, nodes_y, test_function, &residuals_u, &residuals_v);
                cJSON *u_json = cJSON_CreateDoubleArray(residuals_u, RK_ACA->kt);

                cJSON *v_json = cJSON_CreateDoubleArray(residuals_v, RK_ACA->kt);

                cJSON_AddItemToObject(start_json, "residuals_u", u_json);

                cJSON_AddItemToObject(start_json, "residuals_v", v_json);
                rank = RK_ACA->kt;
            }
            else if (strcmp(JSON_STRING(root, "method"), "BACA") == 0)
            {
                RK_ACA = b_aca_rkmatrix(eps_residual, d, dim, rows, cols, nodes_x, nodes_y, test_function, &residuals_u, &residuals_v, &rank_increase);
                cJSON *u_json = cJSON_CreateDoubleArray(residuals_u, RK_ACA->kt);
                cJSON *v_json = cJSON_CreateDoubleArray(residuals_v, RK_ACA->kt);
                cJSON *rank_increase_json = cJSON_CreateDoubleArray(rank_increase, RK_ACA->kt);
                cJSON_AddItemToObject(start_json, "rank_increase", rank_increase_json);
                cJSON_AddItemToObject(start_json, "residuals_u", u_json);
                cJSON_AddItemToObject(start_json, "residuals_v", v_json);
                rank = RK_ACA->kt;
            }
            else if (strcmp(JSON_STRING(root, "method"), "HBACA") == 0)
            {
                int r_h;
                ACAResidualNode *new_residual_node(void);
                ACAResidualNode *root = new_residual_node();
                root->row_start = 0;
                root->col_start = 0;
                h_b_aca_rkmatrix(eps_residual, d, L, root, &U_h, &S_h, &V_h, dim, rows, cols, nodes_x, nodes_y, test_function, &r_h);
                cJSON *tree_json = aca_residual_node_to_json(root);
                cJSON_AddItemToObject(start_json, "h_b_aca_tree", tree_json);
                del_residual_node(root);
                rank = r_h;
            }
            if (strcmp(JSON_STRING(root, "method"), "HBACA") == 0)
            {
                printf("HBACA produced rank %d\n", rank);

                cJSON_AddNumberToObject(start_json, "start", start + 1);

                cJSON_AddNumberToObject(start_json, "rank", rank);
                if (rank >= 3000)
                {
                    printf("Rank %d reached target. No further random starts needed.\n", rank);
                    stop_after_this_run = true;
                }
                cJSON *frob_array = cJSON_CreateArray();

                double *current_frob_norms = malloc(rank * sizeof(double));

                if (!current_frob_norms)
                {
                    printf("ERROR: malloc failed.\n");
                    cJSON_Delete(residual_bench);
                    return 0;
                }

                for (int k = 1; k <= rank; k++)
                {
                    /* truncated U(:,1:k) */
                    pfullmatrix Uk = new_fullmatrix(U_h->rows, k);

                    for (int j = 0; j < k; j++)
                    {
                        for (int i = 0; i < U_h->rows; i++)
                        {
                            Uk->e[j * Uk->rows + i] = U_h->e[j * U_h->rows + i];
                        }
                    }

                    /* truncated S(1:k,1:k) */
                    pfullmatrix Sk = new_zero_fullmatrix(k, k);

                    for (int j = 0; j < k; j++)
                    {
                        for (int i = 0; i < k; i++)
                        {
                            Sk->e[j * Sk->rows + i] = S_h->e[j * S_h->rows + i];
                        }
                    }

                    /* truncated V(1:k,:) */
                    pfullmatrix Vk = new_fullmatrix(k, V_h->cols);

                    for (int j = 0; j < V_h->cols; j++)
                    {
                        for (int i = 0; i < k; i++)
                        {
                            Vk->e[j * Vk->rows + i] = V_h->e[j * V_h->rows + i];
                        }
                    }

                    /* F = Uk*Sk*Vk */
                    pfullmatrix US = new_fullmatrix(Uk->rows, Sk->cols);
                    mul_fullmatrix(Uk, Sk, US);
                    pfullmatrix F = new_fullmatrix(US->rows, Vk->cols);
                    mul_fullmatrix(US, Vk, F);

                    double frob_norm = 0.0;

                    for (int j = 0; j < cols; j++)
                    {
                        for (int i = 0; i < rows; i++)
                        {
                            double xi[dim];
                            double yj[dim];

                            for (int l = 0; l < dim; l++)
                            {
                                xi[l] = nodes_x[l * rows + i];
                                yj[l] = nodes_y[l * cols + j];
                            }

                            double exact = test_function(dim, xi, yj);

                            double approx = F->e[j * F->rows + i];
                            double diff = exact - approx;
                            frob_norm += diff * diff;
                        }
                    }

                    frob_norm = sqrt(frob_norm);
                    current_frob_norms[k - 1] = frob_norm;
                    cJSON_AddItemToArray(frob_array, cJSON_CreateNumber(frob_norm));

                    printf(
                        "HBACA rank %d: "
                        "Frobenius norm = %.6e\n",
                        k, frob_norm);

                    del_fullmatrix(Uk);
                    del_fullmatrix(Sk);
                    del_fullmatrix(Vk);
                    del_fullmatrix(US);
                    del_fullmatrix(F);
                }

                double final_norm = current_frob_norms[rank - 1];

                printf(
                    "HBACA final norm "
                    "for start %d = %.6e\n",
                    start + 1, final_norm);

                cJSON_AddNumberToObject(start_json, "final_norm", final_norm);

                cJSON_AddItemToObject(start_json, "frob_norms", frob_array);

                char key[16];
                snprintf(key, sizeof(key), "%d", start + 1);
                cJSON_AddItemToObject(starts_json, key, start_json);
            }
            else
            {
                printf("%s produced rank %d\n",
                       JSON_STRING(root, "method"),
                       rank);

                if (rank >= 3000)
                {
                    printf("Rank %d reached target. No further random starts needed.\n", rank);
                    stop_after_this_run = true;
                }

                cJSON_AddNumberToObject(start_json, "start", start + 1);
                cJSON_AddNumberToObject(start_json, "rank", rank);

                cJSON *frob_array = cJSON_CreateArray();

                /* determine how many Frobenius norms we will compute */

                int num_frob;

                if (strcmp(JSON_STRING(root, "method"), "ACA") == 0)
                    num_frob = rank;
                else
                {
                    num_frob = 0;

                    while (num_frob < rank &&
                           rank_increase[num_frob] > 0)
                    {
                        num_frob++;
                    }
                }

                double *current_frob_norms =
                    malloc(num_frob * sizeof(double));

                if (!current_frob_norms)
                {
                    printf("ERROR: malloc failed.\n");
                    del_rkmatrix(RK_ACA);
                    cJSON_Delete(residual_bench);
                    return 0;
                }

                RK_ACA->kt = 0;

                if (strcmp(JSON_STRING(root, "method"), "ACA") == 0)
                {
                    RK_ACA->kt = 0;

                    for (int iter = 0; iter < rank; iter++)
                    {
                        RK_ACA->kt = iter + 1;

                        pfullmatrix F =
                            new_zero_fullmatrix(RK_ACA->rows,
                                                RK_ACA->cols);

                        convertrk2_fullmatrix(RK_ACA, F);

                        double frob_norm = 0.0;

                        for (int j = 0; j < RK_ACA->cols; j++)
                        {
                            for (int i = 0; i < RK_ACA->rows; i++)
                            {
                                double xi[dim];
                                double yj[dim];

                                for (int l = 0; l < dim; l++)
                                {
                                    xi[l] = nodes_x[l * rows + i];
                                    yj[l] = nodes_y[l * cols + j];
                                }

                                double exact =
                                    test_function(dim,
                                                  xi,
                                                  yj);

                                double approx =
                                    F->e[j * F->rows + i];

                                double diff =
                                    exact - approx;

                                frob_norm += diff * diff;
                            }
                        }

                        frob_norm = sqrt(frob_norm);

                        current_frob_norms[iter] = frob_norm;

                        cJSON_AddItemToArray(
                            frob_array,
                            cJSON_CreateNumber(frob_norm));

                        printf("ACA rank %d: Frobenius norm = %.6e\n",
                               iter + 1,
                               frob_norm);

                        del_fullmatrix(F);
                    }
                }
                else if (strcmp(JSON_STRING(root, "method"), "BACA") == 0)
                {
                    RK_ACA->kt = 0;

                    int cumulative_rank = 0;

                    for (int iter = 0; iter < num_frob; iter++)
                    {
                        cumulative_rank += (int)rank_increase[iter];

                        RK_ACA->kt = cumulative_rank;

                        pfullmatrix F =
                            new_zero_fullmatrix(RK_ACA->rows,
                                                RK_ACA->cols);

                        convertrk2_fullmatrix(RK_ACA, F);

                        double frob_norm = 0.0;

                        for (int j = 0; j < RK_ACA->cols; j++)
                        {
                            for (int i = 0; i < RK_ACA->rows; i++)
                            {
                                double xi[dim];
                                double yj[dim];

                                for (int l = 0; l < dim; l++)
                                {
                                    xi[l] = nodes_x[l * rows + i];
                                    yj[l] = nodes_y[l * cols + j];
                                }

                                double exact =
                                    test_function(dim,
                                                  xi,
                                                  yj);

                                double approx =
                                    F->e[j * F->rows + i];

                                double diff =
                                    exact - approx;

                                frob_norm += diff * diff;
                            }
                        }

                        frob_norm = sqrt(frob_norm);

                        current_frob_norms[iter] = frob_norm;

                        cJSON_AddItemToArray(
                            frob_array,
                            cJSON_CreateNumber(frob_norm));

                        printf("BACA rank %d: Frobenius norm = %.6e\n",
                               cumulative_rank,
                               frob_norm);

                        del_fullmatrix(F);
                    }
                }

                double final_norm = current_frob_norms[num_frob - 1];

                printf(
                    "ACA final norm "
                    "for start %d = %.6e\n",
                    start + 1, final_norm);

                cJSON_AddNumberToObject(start_json, "final_norm", final_norm);

                cJSON_AddItemToObject(start_json, "frob_norms", frob_array);

                char key[16];
                snprintf(key, sizeof(key), "%d", start + 1);
                cJSON_AddItemToObject(starts_json, key, start_json);

                if (final_norm < best_aca_final_norm)
                {
                    best_aca_final_norm = final_norm;

                    best_aca_rank = rank;

                    best_aca_start = start + 1;

                    if (best_aca_frob_norms)
                    {
                        free(
                            best_aca_frob_norms);
                    }

                    best_aca_frob_norms =
                        malloc(rank * sizeof(double));

                    for (int i = 0; i < num_frob; i++)
                        best_aca_frob_norms[i] = current_frob_norms[i];

                    printf(
                        "--> New best ACA run found.\n");
                }

                free(current_frob_norms);

                del_rkmatrix(RK_ACA);
            }
            if (stop_after_this_run)
            {
                break;
            }
        }
        if (best_aca_frob_norms)
        {
            cJSON *best_json = cJSON_CreateObject();

            cJSON_AddNumberToObject(best_json, "start", best_aca_start);
            cJSON_AddNumberToObject(best_json, "rank", best_aca_rank);
            cJSON_AddNumberToObject(best_json, "final_norm", best_aca_final_norm);

            cJSON *best_frob = cJSON_CreateArray();

            for (int i = 0; i < best_aca_rank; i++)
            {
                cJSON_AddItemToArray(
                    best_frob,
                    cJSON_CreateNumber(best_aca_frob_norms[i]));
            }

            cJSON_AddItemToObject(best_json,
                                  "frob_norms",
                                  best_frob);

            cJSON_AddItemToObject(residual_bench,
                                  "best_run",
                                  best_json);

            free(best_aca_frob_norms);
        }
        cJSON_AddItemToObject(residual_bench,
                              "starts",
                              starts_json);
        cJSON_AddItemToObject(root, "residual_benchmark", residual_bench);
    }

if (hbaca_residual_test)
{
    printf("\n-------------------------------------\n");
    printf("Running HBACA residual test for %s\n",
           JSON_STRING(root, "method"));

    cJSON *hbaca_bench = cJSON_CreateObject();
    cJSON_AddBoolToObject(hbaca_bench, "enabled", true);
    cJSON_AddNumberToObject(hbaca_bench, "epsilon", eps_hbaca);

    cJSON *runs = cJSON_CreateArray();

    for (int iter = 0; iter < h_baca_repeats; iter++)
    {
        printf("\nIteration %d / %d\n", iter + 1, h_baca_repeats);

        cJSON *run = cJSON_CreateObject();

        prkmatrix RK_ACA = NULL;

        pfullmatrix U_h = NULL;
        pfullmatrix S_h = NULL;
        pfullmatrix V_h = NULL;

        double *residuals_u = NULL;
        double *residuals_v = NULL;
        double *rank_increase = NULL;

        int rank = 0;

        if (L > 0)
        {
            printf("HBACA (L=%d, d=%d)\n", L, d);

            ACAResidualNode *root_node = new_residual_node();
            root_node->row_start = 0;
            root_node->col_start = 0;

            h_b_aca_rkmatrix(
                eps_hbaca,
                d,
                L,
                root_node,
                &U_h,
                &S_h,
                &V_h,
                dim,
                rows,
                cols,
                nodes_x,
                nodes_y,
                test_function,
                &rank);

            cJSON_AddItemToObject(
                run,
                "tree",
                aca_residual_node_to_json(root_node));

            del_residual_node(root_node);

            cJSON_AddNumberToObject(run, "rank", rank);

            /* Store singular values returned by HBACA */

            cJSON *sigma_array = cJSON_CreateArray();

            for (int i = 0; i < S_h->rows; i++)
            {
                double sigma = S_h->e[i * S_h->rows + i];

                cJSON_AddItemToArray(
                    sigma_array,
                    cJSON_CreateNumber(sigma));
            }

            cJSON_AddItemToObject(
                run,
                "singular_values",
                sigma_array);

            del_fullmatrix(U_h);
            del_fullmatrix(S_h);
            del_fullmatrix(V_h);
        }
        else if (d > 0)
        {
            printf("BACA (d=%d)\n", d);

            RK_ACA = b_aca_rkmatrix(
                eps_hbaca,
                d,
                dim,
                rows,
                cols,
                nodes_x,
                nodes_y,
                test_function,
                &residuals_u,
                &residuals_v,
                &rank_increase);

            rank = RK_ACA->kt;

            cJSON_AddNumberToObject(run, "rank", rank);

            pfullmatrix A;
            pfullmatrix U;
            pfullmatrix S;
            pfullmatrix V_T;

            A = new_fullmatrix(RK_ACA->rows, RK_ACA->cols);

            convertrk2_fullmatrix(RK_ACA, A);

            int svd_rank = SVD_truncated(
                A,
                1e-14,
                &U,
                &S,
                &V_T);

            cJSON_AddNumberToObject(run, "svd_rank", svd_rank); 
                        /* Store singular values from the truncated SVD */

            cJSON *sigma_array = cJSON_CreateArray();

            for (int i = 0; i < svd_rank; i++)
            {
                double sigma = S->e[i * S->rows + i];

                cJSON_AddItemToArray(
                    sigma_array,
                    cJSON_CreateNumber(sigma));
            }

            cJSON_AddItemToObject(
                run,
                "singular_values",
                sigma_array);

            del_rkmatrix(RK_ACA);

            del_fullmatrix(A);
            del_fullmatrix(U);
            del_fullmatrix(S);
            del_fullmatrix(V_T);

            free(residuals_u);
            free(residuals_v);
            free(rank_increase);
        }
        else
        {
            printf("HBACA residual test requires d > 0.\n");
            cJSON_Delete(run);
            continue;
        }

        cJSON_AddItemToArray(runs, run);
    }

    cJSON_AddItemToObject(hbaca_bench, "runs", runs);
    cJSON_AddItemToObject(root,
                          "hbaca_residual_test",
                          hbaca_bench);
}

    char *json = cJSON_Print(root);
    FILE *fp = fopen(filename, "w");
    if (!fp)
    {
        perror("fopen");
        free(json);
        cJSON_Delete(root);
        return 0;
    }
    fprintf(fp, "%s", json);
    fclose(fp);

    free(json);
    cJSON_Delete(root);
    return 0;
}
// void svd_test(){}
