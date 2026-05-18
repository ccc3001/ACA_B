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

int main()
{
    bool sanity_check = true;
    bool runtime_benchmark = false;
    bool baca_residual_benchmark = false;
    bool h_baca_benchmark = false; // TODO: implement h_baca benchmark
    bool svd_benchmark = false;    // TODO: implement svd benchmark
    /* benchmark settings */
    int iter = 50;
    int repeats = 5;

    /* BACA block sizes */
    int d[] = {2, 4, 8, 16, 32};
    int size_d = sizeof(d) / sizeof(d[0]);

    if (sanity_check)
    {
        printf("=====================================\n");
        printf(" ACA / BACA Sanity Check\n");
        printf("=====================================\n");
        pfullmatrix A_ = gaussian_kernel_matrix(10, 10, 40, 0, 1.0);
        printf("Original matrix A:\n");
        print_fullmatrix(A_);
        prkmatrix RK_ACA = aca_rkmatrix_new(0.00000000000000000001, A_);
        printf("ACA rank: %d\n", RK_ACA->kt);
        pfullmatrix ACA_approx = new_zero_fullmatrix(A_->rows, A_->cols);
        convertrk2_fullmatrix(RK_ACA, ACA_approx);
        printf("ACA approximation:\n");
        print_fullmatrix(ACA_approx);
        prkmatrix RK_BACA = b_aca_rkmatrix_new(0.00000000000000000001, 2, A_);
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

        /* CSV headers */
        fprintf(fbaca, "n,d,elapsed,rank\n");
        fprintf(faca, "n,elapsed,rank\n");

        /* ===================================== */
        /* MAIN BENCHMARK LOOP                   */
        /* ===================================== */

        for (int i = 0; i < iter; i++)
        {
            int n = 1000 + (i + 1) * 100;

            printf("\n=====================================\n");
            printf("Matrix size: %d x %d\n", n, n);
            printf("=====================================\n");

            /* ================================= */
            /* BACA                              */
            /* ================================= */

            for (int j = 0; j < size_d; j++)
            {
                double total_time = 0.0;
                int final_rank = 0;

                printf("\n-------------------------------------\n");
                printf("Running BACA with d = %d\n", d[j]);

                for (int r = 0; r < repeats; r++)
                {
                    printf("repeat %d / %d\n",
                           r + 1,
                           repeats);

                    /* fresh matrix */
                    pfullmatrix A =
                        gaussian_kernel_matrix(
                            n,
                            n,
                            16 * n,
                            0,
                            1.0);

                    double start =
                        omp_get_wtime();

                    prkmatrix RK_BACA =
                        b_aca_rkmatrix_new(
                            0.01,
                            d[j],
                            A);

                    double end =
                        omp_get_wtime();

                    double elapsed =
                        end - start;

                    total_time += elapsed;

                    final_rank =
                        RK_BACA->kt;

                    printf("time : %.6f s\n",
                           elapsed);

                    printf("rank : %d\n",
                           final_rank);

                    del_rkmatrix(RK_BACA);
                    del_fullmatrix(A);
                }

                double avg_time =
                    total_time / repeats;

                printf("\nAverage BACA Results\n");
                printf("n        : %d\n", n);
                printf("d        : %d\n", d[j]);
                printf("avg time : %.6f s\n", avg_time);
                printf("rank     : %d\n", final_rank);

                if (final_rank > 0)
                {
                    printf("T/rank    : %.6e\n",
                           avg_time / final_rank);

                    printf("T/(n*k)   : %.6e\n",
                           avg_time / (n * final_rank));

                    printf("T/(n*k²)  : %.6e\n",
                           avg_time /
                               (n *
                                final_rank *
                                final_rank));
                }

                /* save CSV row */
                fprintf(
                    fbaca,
                    "%d,%d,%.12f,%d\n",
                    n,
                    d[j],
                    avg_time,
                    final_rank);
            }

            /* ================================= */
            /* ACA                               */
            /* ================================= */

            double total_time_ACA = 0.0;
            int final_rank_ACA = 0;

            printf("\n-------------------------------------\n");
            printf("Running ACA\n");

            for (int r = 0; r < repeats; r++)
            {
                printf("repeat %d / %d\n",
                       r + 1,
                       repeats);

                pfullmatrix A =
                    gaussian_kernel_matrix(
                        n,
                        n,
                        16 * n,
                        0,
                        1.0);

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

                total_time_ACA += elapsed;

                final_rank_ACA =
                    RK_ACA->kt;

                printf("time : %.6f s\n",
                       elapsed);

                printf("rank : %d\n",
                       final_rank_ACA);

                del_rkmatrix(RK_ACA);
                del_fullmatrix(A);
            }

            double avg_time_ACA =
                total_time_ACA / repeats;

            printf("\nAverage ACA Results\n");
            printf("n        : %d\n", n);
            printf("avg time : %.6f s\n",
                   avg_time_ACA);

            printf("rank     : %d\n",
                   final_rank_ACA);

            if (final_rank_ACA > 0)
            {
                printf("T/rank    : %.6e\n",
                       avg_time_ACA /
                           final_rank_ACA);

                printf("T/(n*k)   : %.6e\n",
                       avg_time_ACA /
                           (n * final_rank_ACA));

                printf("T/(n*k²)  : %.6e\n",
                       avg_time_ACA /
                           (n *
                            final_rank_ACA *
                            final_rank_ACA));
            }

            /* save CSV row */
            fprintf(
                faca,
                "%d,%.12f,%d\n",
                n,
                avg_time_ACA,
                final_rank_ACA);
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

    if (baca_residual_benchmark)
    {
        printf("=====================================\n");
        printf(" BACA Frobenius Norm over time\n");
        printf("=====================================\n");
        FILE *frob = fopen("frobenius_norms.csv", "w");
        if (!frob)
        {
            printf("ERROR: could not open output files.\n");
            return 1;
        }
        fprintf(frob, "d,rank,frobenius_norm\n");
        pfullmatrix A =
            gaussian_kernel_matrix(
                5000,
                5000,
                2 * 5000,
                0,
                1.0);

        prkmatrix RK_ACA = aca_rkmatrix_new(
            0.001,
            A);

        int kt = RK_ACA->kt;
        double *frob_norms = malloc(kt * sizeof(double));
        RK_ACA->kt = 0;
        for (int i = 1; i <= kt; i++)
        {
            RK_ACA->kt = i;
            pfullmatrix F = new_zero_fullmatrix(A->rows, A->cols);
            convertrk2_fullmatrix(RK_ACA, F);
            double frob_norm = 0.0;
            for (int k = 0; k < A->rows * A->cols; k++)
            {
                double diff = A->e[k] - F->e[k];
                frob_norm += sqrt(diff * diff);
            }
            printf("ACA rank %d: Frobenius norm of (A - R) = %.6e\n", i, frob_norm);
            frob_norms[i - 1] = frob_norm;
            del_fullmatrix(F);
        }
        for (int i = 0; i < kt; i++)
        {
            fprintf(frob, "%d,%d,%.6e\n", 1, i + 1, frob_norms[i]);
        }

        for (int j = 0; j < size_d; j++)
        {

            prkmatrix RK_BACA = b_aca_rkmatrix_new(
                0.001,
                d[j],
                A);

            printf("d = %d, rank = %d\n", d[j], RK_BACA->kt);
            int kt = RK_BACA->kt;
            double *frob_norms = malloc(kt * sizeof(double));
            RK_BACA->kt = 0;
            for (int i = 1; i <= kt; i++)
            {
                RK_BACA->kt = i;
                pfullmatrix F = new_zero_fullmatrix(A->rows, A->cols);
                convertrk2_fullmatrix(RK_BACA, F);
                double frob_norm = 0.0;
                for (int k = 0; k < A->rows * A->cols; k++)
                {
                    double diff = A->e[k] - F->e[k];
                    frob_norm += sqrt(diff * diff);
                }
                printf("rank %d: Frobenius norm of (A - R) = %.6e\n", i, frob_norm);
                frob_norms[i - 1] = frob_norm;
                del_fullmatrix(F);
                del_rkmatrix(RK_BACA);
            }
            for (int i = 0; i < kt; i++)
            {
                fprintf(frob, "%d,%d,%.6e\n", d[j], i + 1, frob_norms[i]);
            }
        }
        fclose(frob);
        del_doubles(frob_norms);
        del_fullmatrix(A);
        del_rkmatrix(RK_ACA);
    }

    return 0;
}