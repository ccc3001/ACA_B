#include "new_aca.h"
#include "aca_b.h"
#include "fullmatrix.h"
#include "rkmatrix.h"
#include "basic.h"
#include "kernel_functions.h"
#include "svd.h"

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>

int main()
{
    /* BACA block sizes */
    int d[] = {2, 4, 8, 16, 32};
    int size_d = sizeof(d) / sizeof(d[0]);

    /*H_BACA depth L*/
    int L[] = {4, 5};
    int size_L = sizeof(L) / sizeof(L[0]);

    bool sanity_check = false;
    bool runtime_benchmark = true;
    /* runtime benchmark settings */
    int iter = 20;
    int repeats = 4;

    /* residual benchmark settings */
    int n_residual = 1000;
    double h_residual = 1.0;
    double eps_residual = 0.00001;

    bool baca_residual_benchmark = false;

    bool h_baca_benchmark = false; // TODO: implement h_baca benchmark
    bool svd_benchmark = false;    // TODO: implement svd benchmark

    if (sanity_check)
    {
        printf("=====================================\n");
        printf(" ACA / BACA Sanity Check\n");
        printf("=====================================\n");
        pfullmatrix A_ = gaussian_kernel_matrix(10, 10, 40, 0, 1.0);
        printf("Original matrix A:\n");
        print_fullmatrix(A_);
        double *residuals;
        prkmatrix RK_ACA = aca_rkmatrix_new(0.00000000000000000001, A_);
        printf("ACA rank: %d\n", RK_ACA->kt);
        pfullmatrix ACA_approx = new_zero_fullmatrix(A_->rows, A_->cols);
        convertrk2_fullmatrix(RK_ACA, ACA_approx);
        printf("ACA approximation:\n");
        print_fullmatrix(ACA_approx);
        prkmatrix RK_BACA = b_aca_rkmatrix_new(0.00000000000000000001, 2, A_, &residuals);
        free(residuals);
        printf("BACA rank: %d\n", RK_BACA->kt);
        pfullmatrix BACA_approx = new_zero_fullmatrix(A_->rows, A_->cols);
        convertrk2_fullmatrix(RK_BACA, BACA_approx);
        printf("BACA approximation:\n");
        print_fullmatrix(BACA_approx);

        del_fullmatrix(ACA_approx);
        del_fullmatrix(BACA_approx);
        del_rkmatrix(RK_ACA);
        del_rkmatrix(RK_BACA);

        printf("=====================================\n");
        printf(" SVD Sanity Check\n");
        printf("=====================================\n");
        pfullmatrix U, S, V_T;
        SVD(A_, &U, &S, &V_T);
        printf("U:\n");
        print_fullmatrix(U);
        printf("S:\n");
        print_fullmatrix(S);
        printf("V^T:\n");
        print_fullmatrix(V_T);
        pfullmatrix US = new_zero_fullmatrix(U->rows, S->cols);
        mul_fullmatrix(U, S, US);
        pfullmatrix SVD_approx = new_zero_fullmatrix(A_->rows, A_->cols);
        mul_fullmatrix(US, V_T, SVD_approx);
        printf("SVD approximation:\n");
        print_fullmatrix(SVD_approx);
        del_fullmatrix(US);
        del_fullmatrix(SVD_approx);
        del_fullmatrix(U);
        del_fullmatrix(S);
        del_fullmatrix(V_T);

        printf("=====================================\n");
        printf(" Submatrix Sanity Check\n");
        printf("=====================================\n");
        pfullmatrix sub11 = new_submatrix(A_, 0, 0, A_->rows / 2, A_->cols / 2);
        printf("Submatrix (1,1):\n");
        print_fullmatrix(sub11);
        del_fullmatrix(sub11);
        pfullmatrix sub12 = new_submatrix(A_, 0, A_->cols / 2, A_->rows / 2, A_->cols / 2);
        printf("Submatrix (1,2):\n");
        print_fullmatrix(sub12);
        del_fullmatrix(sub12);
        pfullmatrix sub21 = new_submatrix(A_, A_->rows / 2, 0, A_->rows / 2, A_->cols / 2);
        printf("Submatrix (2,1):\n");
        print_fullmatrix(sub21);
        del_fullmatrix(sub21);
        pfullmatrix sub22 = new_submatrix(A_, A_->rows / 2, A_->cols / 2, A_->rows / 2, A_->cols / 2);
        printf("Submatrix (2,2):\n");
        print_fullmatrix(sub22);
        del_fullmatrix(sub22);

        printf("=====================================\n");
        printf(" H_BACA Sanity Check\n");
        printf("=====================================\n");

        pfullmatrix U_hbaca, S_hbaca, V_hbaca;
        int r_hbaca;
        h_b_aca_rkmatrix_new(0.00000000000000000001, 2, 2, A_, &U_hbaca, &S_hbaca, &V_hbaca, &r_hbaca);
        printf("H_BACA rank: %d\n", r_hbaca);
        printf("U_hbaca:\n");
        print_fullmatrix(U_hbaca);
        printf("S_hbaca:\n");
        print_fullmatrix(S_hbaca);
        printf("V_hbaca:\n");
        print_fullmatrix(V_hbaca);
        pfullmatrix US_hbaca = new_zero_fullmatrix(U_hbaca->rows, S_hbaca->cols);
        mul_fullmatrix(U_hbaca, S_hbaca, US_hbaca);
        printf("A:  \n");
        print_fullmatrix(A_);
        pfullmatrix H_BACA_approx = new_zero_fullmatrix(A_->rows, A_->cols);

        mul_fullmatrix(US_hbaca, V_hbaca, H_BACA_approx);

        printf("H_BACA approximation:\n");
        print_fullmatrix(H_BACA_approx);
        printf("Difference (A - H_BACA approximation):\n");
        print_matrix_difference(A_, H_BACA_approx);
        del_fullmatrix(US_hbaca);
        del_fullmatrix(H_BACA_approx);
        del_fullmatrix(U_hbaca);
        del_fullmatrix(S_hbaca);
        del_fullmatrix(V_hbaca);

        del_fullmatrix(A_);
    }

    if (runtime_benchmark)
    {
        printf("=====================================\n");
        printf(" ACA / BACA Benchmark\n");
        printf("=====================================\n");

        /* ---------------------------------- */
        /* open output files                  */
        /* ---------------------------------- */

        FILE *fbaca = fopen("baca_results.csv", "w");
        FILE *faca = fopen("aca_results.csv", "w");

        if (!fbaca || !faca)
        {
            printf("ERROR: could not open output files.\n");
            return 1;
        }

        /*
        Save EVERY repeat individually.
        This is MUCH better statistically.
        */

        fprintf(
            fbaca,
            "n,d,repeat,elapsed,rank,T_over_rank,T_over_nk,T_over_nk2\n");

        fprintf(
            faca,
            "n,repeat,elapsed,rank,T_over_rank,T_over_nk,T_over_nk2\n");

        /* ===================================== */
        /* MAIN BENCHMARK LOOP                   */
        /* ===================================== */

        for (int i = 0; i < iter; i++)
        {
            int n = 1000 + (i + 1) * 100;

            printf("\n=====================================\n");
            printf("Matrix size: %d x %d\n", n, n);
            printf("=====================================\n");

            /*
            IMPORTANT:
            Generate ONE matrix per n.
            This isolates algorithmic randomness
            from matrix randomness.
            */

            pfullmatrix A =
                gaussian_kernel_matrix(
                    n,
                    n,
                    16 * n,
                    0,
                    1.0);

            /* ================================= */
            /* BACA                              */
            /* ================================= */

            for (int j = 0; j < size_d; j++)
            {
                double total_time = 0.0;
                double total_rank = 0.0;

                double total_T_over_rank = 0.0;
                double total_T_over_nk = 0.0;
                double total_T_over_nk2 = 0.0;

                printf("\n-------------------------------------\n");
                printf("Running BACA with d = %d\n", d[j]);

                /*
                 * Parallel repeats
                 */

#pragma omp parallel for reduction(+ : total_time, total_rank, total_T_over_rank, total_T_over_nk, total_T_over_nk2)
                for (int r = 0; r < repeats; r++)
                {
                    /*
                     * Per-thread deterministic seed
                     */

                    unsigned int seed =
                        12345u + 1000u * omp_get_thread_num() + (unsigned int)r;

                    /*
                     * Seed RNG for this thread/run
                     */

                    srand(seed);

#pragma omp critical
                    {
                        printf("repeat %d / %d\n",
                               r + 1,
                               repeats);
                    }

                    double start =
                        omp_get_wtime();
                    double *residuals;
                    prkmatrix RK_BACA =
                        b_aca_rkmatrix_new(
                            0.01,
                            d[j],
                            A, &residuals);
                    free(residuals);

                    double end =
                        omp_get_wtime();

                    double elapsed =
                        end - start;

                    int rank =
                        RK_BACA->kt;

                    double T_over_rank = 0.0;
                    double T_over_nk = 0.0;
                    double T_over_nk2 = 0.0;

                    if (rank > 0)
                    {
                        T_over_rank =
                            elapsed / rank;

                        T_over_nk =
                            elapsed /
                            (n * rank);

                        T_over_nk2 =
                            elapsed /
                            (n * rank * rank);
                    }

                    /*
                     * Reduction variables
                     */

                    total_time += elapsed;
                    total_rank += rank;

                    total_T_over_rank += T_over_rank;
                    total_T_over_nk += T_over_nk;
                    total_T_over_nk2 += T_over_nk2;

                    /*
                     * Console output
                     */

#pragma omp critical
                    {
                        printf("time : %.6f s\n",
                               elapsed);

                        printf("rank : %d\n",
                               rank);

                        printf("T/rank   : %.6e\n",
                               T_over_rank);

                        printf("T/(n*k)  : %.6e\n",
                               T_over_nk);

                        printf("T/(n*k²) : %.6e\n",
                               T_over_nk2);
                    }

                    /*
                     * Save EVERY run individually
                     */

#pragma omp critical
                    {
                        fprintf(
                            fbaca,
                            "%d,%d,%d,%.12f,%d,%.12e,%.12e,%.12e\n",
                            n,
                            d[j],
                            r,
                            elapsed,
                            rank,
                            T_over_rank,
                            T_over_nk,
                            T_over_nk2);
                    }

                    del_rkmatrix(RK_BACA);
                }

                /*
                 * Aggregate statistics
                 */

                double avg_time =
                    total_time / repeats;

                double avg_rank =
                    total_rank / repeats;

                double avg_T_over_rank =
                    total_T_over_rank / repeats;

                double avg_T_over_nk =
                    total_T_over_nk / repeats;

                double avg_T_over_nk2 =
                    total_T_over_nk2 / repeats;

                printf("\nAverage BACA Results\n");

                printf("n            : %d\n", n);
                printf("d            : %d\n", d[j]);

                printf("avg time     : %.6f s\n",
                       avg_time);

                printf("avg rank     : %.6f\n",
                       avg_rank);

                printf("avg T/rank   : %.6e\n",
                       avg_T_over_rank);

                printf("avg T/(n*k)  : %.6e\n",
                       avg_T_over_nk);

                printf("avg T/(n*k²) : %.6e\n",
                       avg_T_over_nk2);
            }

            /* ================================= */
            /* ACA                               */
            /* ================================= */

            {
                double total_time = 0.0;
                double total_rank = 0.0;

                double total_T_over_rank = 0.0;
                double total_T_over_nk = 0.0;
                double total_T_over_nk2 = 0.0;

                printf("\n-------------------------------------\n");
                printf("Running ACA\n");

                for (int r = 0; r < repeats; r++)
                {
                    printf("repeat %d / %d\n",
                           r + 1,
                           repeats);

                    double start =
                        omp_get_wtime();

                    prkmatrix RK_ACA =
                        aca_rkmatrix_new(
                            0.01,
                            A);

                    double end =
                        omp_get_wtime();

                    double elapsed =
                        end - start;

                    int rank =
                        RK_ACA->kt;

                    double T_over_rank = 0.0;
                    double T_over_nk = 0.0;
                    double T_over_nk2 = 0.0;

                    if (rank > 0)
                    {
                        T_over_rank =
                            elapsed / rank;

                        T_over_nk =
                            elapsed /
                            (n * rank);

                        T_over_nk2 =
                            elapsed /
                            (n * rank * rank);
                    }

                    total_time += elapsed;
                    total_rank += rank;

                    total_T_over_rank += T_over_rank;
                    total_T_over_nk += T_over_nk;
                    total_T_over_nk2 += T_over_nk2;

                    printf("time : %.6f s\n",
                           elapsed);

                    printf("rank : %d\n",
                           rank);

                    printf("T/rank   : %.6e\n",
                           T_over_rank);

                    printf("T/(n*k)  : %.6e\n",
                           T_over_nk);

                    printf("T/(n*k²) : %.6e\n",
                           T_over_nk2);

                    /*
                    Save EVERY run individually
                    */

                    fprintf(
                        faca,
                        "%d,%d,%.12f,%d,%.12e,%.12e,%.12e\n",
                        n,
                        r,
                        elapsed,
                        rank,
                        T_over_rank,
                        T_over_nk,
                        T_over_nk2);

                    del_rkmatrix(RK_ACA);
                }

                /*
                Aggregate statistics
                */

                double avg_time =
                    total_time / repeats;

                double avg_rank =
                    total_rank / repeats;

                double avg_T_over_rank =
                    total_T_over_rank / repeats;

                double avg_T_over_nk =
                    total_T_over_nk / repeats;

                double avg_T_over_nk2 =
                    total_T_over_nk2 / repeats;

                printf("\nAverage ACA Results\n");

                printf("n            : %d\n", n);

                printf("avg time     : %.6f s\n",
                       avg_time);

                printf("avg rank     : %.6f\n",
                       avg_rank);

                printf("avg T/rank   : %.6e\n",
                       avg_T_over_rank);

                printf("avg T/(n*k)  : %.6e\n",
                       avg_T_over_nk);

                printf("avg T/(n*k²) : %.6e\n",
                       avg_T_over_nk2);
            }

            /*
            cleanup matrix
            */

            del_fullmatrix(A);
        }

        /* ===================================== */
        /* CLOSE FILES                           */
        /* ===================================== */

        fclose(fbaca);
        fclose(faca);

        printf("\n=====================================\n");
        printf("Benchmark complete.\n");
        printf("Saved:\n");
        printf("  - baca_results.csv\n");
        printf("  - aca_results.csv\n");
        printf("=====================================\n");
    }

    pfullmatrix A =
        gaussian_kernel_matrix(
            n_residual,
            n_residual,
            2 * n_residual,
            0,
            h_residual);

    if (baca_residual_benchmark)
    {
        printf("=====================================\n");
        printf(" BACA Frobenius Norm over time\n");
        printf("=====================================\n");

        FILE *frob = fopen("frobenius_norms.csv", "w");

        if (!frob)
        {
            printf("ERROR: could not open output file.\n");
            return 1;
        }

        fprintf(frob, "method,d,start,rank,frobenius_norm\n");

        /*
         * Number of random restarts
         */
        int n_starts = 10;

        /*
         * ============================================================
         * ACA MULTI-START BENCHMARK
         * ============================================================
         */

        printf("\n=====================================\n");
        printf(" ACA MULTI-START BENCHMARK\n");
        printf("=====================================\n");

        double best_aca_final_norm = HUGE_VAL;

        double *best_aca_frob_norms = NULL;

        int best_aca_rank = 0;
        int best_aca_start = -1;

        for (int start = 0; start < n_starts; start++)
        {
            printf(
                "\nACA random start %d / %d\n",
                start + 1,
                n_starts);

            /*
             * Seed RNG for different starts
             */
            srand(time(NULL) + start);

            prkmatrix RK_ACA =
                aca_rkmatrix_new(
                    eps_residual,
                    A);

            int kt_aca = RK_ACA->kt;

            printf(
                "ACA produced rank %d\n",
                kt_aca);

            double *current_frob_norms =
                malloc(kt_aca * sizeof(double));

            if (!current_frob_norms)
            {
                printf("ERROR: malloc failed.\n");

                del_rkmatrix(RK_ACA);
                fclose(frob);

                return 1;
            }

            RK_ACA->kt = 0;

            for (int i = 1; i <= kt_aca; i++)
            {
                RK_ACA->kt = i;

                pfullmatrix F =
                    new_zero_fullmatrix(
                        A->rows,
                        A->cols);

                convertrk2_fullmatrix(
                    RK_ACA,
                    F);

                double frob_norm = 0.0;

                for (int k = 0; k < A->rows * A->cols; k++)
                {
                    double diff =
                        A->e[k] - F->e[k];

                    frob_norm += diff * diff;
                }

                frob_norm = sqrt(frob_norm);

                current_frob_norms[i - 1] =
                    frob_norm;

                printf(
                    "ACA rank %d: Frobenius norm = %.6e\n",
                    i,
                    frob_norm);

                del_fullmatrix(F);
            }

            double final_norm =
                current_frob_norms[kt_aca - 1];

            printf(
                "ACA final norm for start %d = %.6e\n",
                start + 1,
                final_norm);

            /*
             * Keep best ACA run
             */
            if (final_norm < best_aca_final_norm)
            {
                best_aca_final_norm = final_norm;

                best_aca_rank = kt_aca;

                best_aca_start = start + 1;

                if (best_aca_frob_norms)
                {
                    free(best_aca_frob_norms);
                }

                best_aca_frob_norms =
                    malloc(kt_aca * sizeof(double));

                for (int i = 0; i < kt_aca; i++)
                {
                    best_aca_frob_norms[i] =
                        current_frob_norms[i];
                }

                printf(
                    "--> New best ACA run found.\n");
            }

            free(current_frob_norms);

            del_rkmatrix(RK_ACA);
        }

        /*
         * Write best ACA trajectory
         */
        printf(
            "\nBest ACA run: start %d, final norm %.6e\n",
            best_aca_start,
            best_aca_final_norm);

        for (int i = 0; i < best_aca_rank; i++)
        {
            fprintf(
                frob,
                "ACA,%d,%d,%d,%.6e\n",
                1,
                best_aca_start,
                i + 1,
                best_aca_frob_norms[i]);
        }

        free(best_aca_frob_norms);

        /*
         * ============================================================
         * BACA MULTI-START BENCHMARK
         * ============================================================
         */

        printf("\n=====================================\n");
        printf(" BACA MULTI-START BENCHMARK\n");
        printf("=====================================\n");

        for (int j = 0; j < size_d; j++)
        {
            printf(
                "\nStarting calculation for d = %d\n",
                d[j]);

            double best_baca_final_norm = HUGE_VAL;

            double *best_baca_frob_norms = NULL;

            int best_baca_rank = 0;

            int best_baca_start = -1;

            for (int start = 0; start < n_starts; start++)
            {
                printf(
                    "\nBACA random start %d / %d\n",
                    start + 1,
                    n_starts);

                /*
                 * Seed RNG differently each run
                 */
                srand(time(NULL) + start);
                double *residuals;
                prkmatrix RK_BACA =
                    b_aca_rkmatrix_new(
                        eps_residual,
                        d[j],
                        A, &residuals);
                free(residuals);
                int kt_baca = RK_BACA->kt;

                printf(
                    "BACA produced rank %d\n",
                    kt_baca);

                double *current_frob_norms =
                    malloc(kt_baca * sizeof(double));

                if (!current_frob_norms)
                {
                    printf("ERROR: malloc failed.\n");

                    del_rkmatrix(RK_BACA);

                    if (best_baca_frob_norms)
                    {
                        free(best_baca_frob_norms);
                    }

                    fclose(frob);

                    return 1;
                }

                RK_BACA->kt = 0;

                for (int i = 1; i <= kt_baca; i++)
                {
                    RK_BACA->kt = i;

                    pfullmatrix F =
                        new_zero_fullmatrix(
                            A->rows,
                            A->cols);

                    convertrk2_fullmatrix(
                        RK_BACA,
                        F);

                    double frob_norm = 0.0;

                    for (int k = 0; k < A->rows * A->cols; k++)
                    {
                        double diff =
                            A->e[k] - F->e[k];

                        frob_norm += diff * diff;
                    }

                    frob_norm = sqrt(frob_norm);

                    current_frob_norms[i - 1] =
                        frob_norm;

                    printf(
                        "BACA rank %d: Frobenius norm = %.6e\n",
                        i,
                        frob_norm);

                    del_fullmatrix(F);
                }

                double final_norm =
                    current_frob_norms[kt_baca - 1];

                printf(
                    "BACA final norm for start %d = %.6e\n",
                    start + 1,
                    final_norm);

                /*
                 * Keep best BACA run
                 */
                if (final_norm < best_baca_final_norm)
                {
                    best_baca_final_norm = final_norm;

                    best_baca_rank = kt_baca;

                    best_baca_start = start + 1;

                    if (best_baca_frob_norms)
                    {
                        free(best_baca_frob_norms);
                    }

                    best_baca_frob_norms =
                        malloc(kt_baca * sizeof(double));

                    for (int i = 0; i < kt_baca; i++)
                    {
                        best_baca_frob_norms[i] =
                            current_frob_norms[i];
                    }

                    printf(
                        "--> New best BACA run found.\n");
                }

                free(current_frob_norms);

                del_rkmatrix(RK_BACA);
            }

            /*
             * Write best BACA trajectory
             */
            printf(
                "\nBest BACA run for d = %d:\n",
                d[j]);

            printf(
                "Best start = %d\n",
                best_baca_start);

            printf(
                "Best final Frobenius norm = %.6e\n",
                best_baca_final_norm);

            for (int i = 0; i < best_baca_rank; i++)
            {
                fprintf(
                    frob,
                    "BACA,%d,%d,%d,%.6e\n",
                    d[j],
                    best_baca_start,
                    i + 1,
                    best_baca_frob_norms[i]);
            }

            free(best_baca_frob_norms);
        }

        fclose(frob);
    }
    if (svd_benchmark)
    {
        FILE *frob = fopen("frobenius_norms_svd.csv", "w");

        if (!frob)
        {
            printf("ERROR: could not open output files.\n");
            return 1;
        }

        fprintf(frob, "rank,frobenius_norm\n");

        pfullmatrix U, V_T, S;

        SVD(A, &U, &S, &V_T);

        int r = S->rows;

        /* -----------------------------------------
        Total squared Frobenius norm
        ----------------------------------------- */

        double tail_sum = 0.0;

        for (int i = 0; i < r; i++)
        {
            double sigma = get_fullmatrix_value(S, i, i);

            tail_sum += sigma * sigma;
        }

        /* -----------------------------------------
        Rank-k truncation errors
        ----------------------------------------- */

        for (int k = 0; k < r; k++)
        {
            double sigma = get_fullmatrix_value(S, k, k);

            tail_sum -= sigma * sigma;

            double frob_norm = sqrt(fmax(tail_sum, 0.0));

            printf(
                "SVD rank %d: Frobenius norm of (A - A_k) = %.6e\n",
                k + 1,
                frob_norm);

            fprintf(
                frob,
                "%d,%.6e\n",
                k + 1,
                frob_norm);
        }

        fclose(frob);

        del_fullmatrix(U);
        del_fullmatrix(S);
        del_fullmatrix(V_T);
    }
    if (h_baca_benchmark)
    {
        for (int j = 0; j < size_d; j++)
        {
            for (int i = 0; i < size_L; i++)
            {
                pfullmatrix U, V, S;
                int r;
                h_b_aca_rkmatrix_new(eps_residual, d[j], L[i], A, &U, &S, &V, &r);
            }
        }

        printf("work in progress");
    }
    del_fullmatrix(A);

    return 0;
}