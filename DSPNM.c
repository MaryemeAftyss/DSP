#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// Fonction pour générer un bruit gaussien
double bruit_gaussien(double moyenne, double ecart_type) {
    double u1 = rand() / (double)RAND_MAX;
    double u2 = rand() / (double)RAND_MAX;
    double z0 = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
    return moyenne + z0 * ecart_type;
}

// Fonction pour la démodulation (trouve le symbole le plus proche)
int demodulation(double x, double y, int N) {
    int i, j;
    double distance_min = __DBL_MAX__;
    int symbole_demodule = -1;

    for (i = -(N-1); i <= (N-1); i += 2) {
        for (j = -(N-1); j <= (N-1); j += 2) {
            double distance = sqrt(pow(x - i, 2) + pow(y - j, 2));
            if (distance < distance_min) {
                distance_min = distance;
                symbole_demodule = i * N + j;
            }
        }
    }

    return symbole_demodule;
}

int main() {
    int m;
    printf("Entrez le nombre de bits par symbole (m) : ");
    scanf("%d", &m);

    if (m <= 0) {
        printf("Le nombre de bits par symbole doit être positif.\n");
        return 1;
    }

    int M = (int)pow(2, m); // Nombre total de symboles

    printf("Choisissez le type de bruit :\n");
    printf("1. Bruit additif gaussien (AWGN)\n");
    printf("2. Bruit uniforme\n");
    printf("Votre choix : ");
    int bruit_type;
    scanf("%d", &bruit_type);

    double bruit_intensite;
    printf("Entrez l'intensité du bruit (par exemple 0.5) : ");
    scanf("%lf", &bruit_intensite);

    srand(time(NULL));

    // Variables pour le calcul du BER
    int erreurs = 0;
    int total_bits = 0;

    // Cas où m est pair
    if (m % 2 == 0) {
        int k = m / 2;
        int N = (int)pow(2, k); // Nombre de lignes et colonnes

        printf("Pour %d-QAM (m = %d bits par symbole, pair) :\n", M, m);
        printf("On aura %d lignes et %d colonnes\n\n", N, N);

        // Variables pour le calcul du SNR
        double total_signal_power = 0.0;
        double total_noise_power = 0.0;
        int symbol_count = 0;

        // Affichage des points de la constellation
        int symbol_index = 0;
        for (int row = -(N - 1); row <= (N - 1); row += 2) {
            for (int col = -(N - 1); col <= (N - 1); col += 2) {
                double noisy_row = row;
                double noisy_col = col;

                // Ajouter du bruit aux points de la constellation
                if (bruit_type == 1) {
                    noisy_row += bruit_gaussien(0, bruit_intensite);
                    noisy_col += bruit_gaussien(0, bruit_intensite);
                } else if (bruit_type == 2) {
                    noisy_row += ((rand() / (double)RAND_MAX) * 2.0 - 1.0) * bruit_intensite;
                    noisy_col += ((rand() / (double)RAND_MAX) * 2.0 - 1.0) * bruit_intensite;
                }

                // Calcul de la puissance du signal
                total_signal_power += pow(row, 2) + pow(col, 2);

                // Calcul de la puissance du bruit
                total_noise_power += pow(noisy_row - row, 2) + pow(noisy_col - col, 2);

                // Démodulation : récupération du symbole le plus proche
                int demod_symbol = demodulation(noisy_row, noisy_col, N);

                // Comparer les bits transmis et démoulus
                if (symbol_index != demod_symbol) {
                    erreurs++; // Incrémenter les erreurs si les symboles sont différents
                }

                total_bits++; // Compter le nombre total de bits transmis

                // Affichage des points de la constellation et des erreurs
                printf("Point original (%d, %d) -> Avec bruit (%.2f, %.2f): Symbole %d (démodulé %d)\n", row, col, noisy_row, noisy_col, symbol_index, demod_symbol);
                symbol_index++;
            }
        }

        // Calcul du taux d'erreur binaire (BER)
        double ber = (double)erreurs / total_bits;
        printf("\nLe taux d'erreur binaire (BER) est : %.5f\n", ber);

        // Calcul et affichage du SNR
        if (total_noise_power > 0) {
            double snr_db = 10 * log10(total_signal_power / total_noise_power);
            printf("\nSNR moyen = %.2f dB\n", snr_db);
        } else {
            printf("\nLe bruit est nul, impossible de calculer le SNR.\n");
        }

        // Bilan énergétique
        printf("\nBilan énergétique :\n");
        for (int row = -(N - 1); row <= (N - 1); row += 2) {
            for (int col = -(N - 1); col <= (N - 1); col += 2) {
                double energy = pow(row, 2) + pow(col, 2);
                printf("Point (%d, %d): Energie = %.2f\n", row, col, energy);
            }
        }

        // Bilan de phase
        printf("\nBilan de phase :\n");
        for (int row = -(N - 1); row <= (N - 1); row += 2) {
            for (int col = -(N - 1); col <= (N - 1); col += 2) {
                double angle = atan2(col, row);
                printf("Point (%d, %d): Phase = %.2f rad, %.2f°\n", row, col, angle, angle * (180 / M_PI));
            }
        }
    } else {
        // Cas où m est impair
        int N = (int)ceil(sqrt(M)); // Nombre de lignes et colonnes nécessaires
        int symbols_to_remove = N * N - M; // Symboles à exclure

        printf("Pour %d-QAM (m = %d bits par symbole, impair) :\n", M, m);
        printf("On aura %d lignes et %d colonnes en éliminant %d symboles.\n\n", N, N, symbols_to_remove);

        // Variables pour le calcul du SNR
        double total_signal_power = 0.0;
        double total_noise_power = 0.0;
        int symbol_count = 0;

        // Affichage des points de la constellation
        int symbol_index = 0;
        for (int row = -(N - 1); row <= (N - 1); row += 2) {
            for (int col = -(N - 1); col <= (N - 1); col += 2) {
                if (symbol_index < M) {
                    double noisy_row = row;
                    double noisy_col = col;

                    // Ajouter du bruit aux points de la constellation
                    if (bruit_type == 1) {
                        noisy_row += bruit_gaussien(0, bruit_intensite);
                        noisy_col += bruit_gaussien(0, bruit_intensite);
                    } else if (bruit_type == 2) {
                        noisy_row += ((rand() / (double)RAND_MAX) * 2.0 - 1.0) * bruit_intensite;
                        noisy_col += ((rand() / (double)RAND_MAX) * 2.0 - 1.0) * bruit_intensite;
                    }

                    // Calcul de la puissance du signal
                    total_signal_power += pow(row, 2) + pow(col, 2);

                    // Calcul de la puissance du bruit
                    total_noise_power += pow(noisy_row - row, 2) + pow(noisy_col - col, 2);

                    // Démodulation : récupération du symbole le plus proche
                    int demod_symbol = demodulation(noisy_row, noisy_col, N);

                    // Comparer les bits transmis et démoulus
                    if (symbol_index != demod_symbol) {
                        erreurs++; // Incrémenter les erreurs si les symboles sont différents
                    }

                    total_bits++; // Compter le nombre total de bits transmis

                    // Affichage des points de la constellation et des erreurs
                    printf("Point original (%d, %d) -> Avec bruit (%.2f, %.2f): Symbole %d (démodulé %d)\n", row, col, noisy_row, noisy_col, symbol_index, demod_symbol);
                    symbol_index++;
                }
            }
        }

        // Calcul du taux d'erreur binaire (BER)
        double ber = (double)erreurs / total_bits;
        printf("\nLe taux d'erreur binaire (BER) est : %.5f\n", ber);

        // Calcul et affichage du SNR
        if (total_noise_power > 0) {
            double snr_db = 10 * log10(total_signal_power / total_noise_power);
            printf("\nSNR moyen = %.2f dB\n", snr_db);
        } else {
            printf("\nLe bruit est nul, impossible de calculer le SNR.\n");
        }

        // Bilan énergétique
        printf("\nBilan énergétique :\n");
        symbol_index = 0;
        for (int row = -(N - 1); row <= (N - 1); row += 2) {
            for (int col = -(N - 1); col <= (N - 1); col += 2) {
                if (symbol_index < M) {
                    double energy = pow(row, 2) + pow(col, 2);
                    printf("Point (%d, %d): Energie = %.2f\n", row, col, energy);
                    symbol_index++;
                }
            }
        }

        // Bilan de phase
        printf("\nBilan de phase :\n");
        symbol_index = 0;
        for (int row = -(N - 1); row <= (N - 1); row += 2) {
            for (int col = -(N - 1); col <= (N - 1); col += 2) {
                if (symbol_index < M) {
                    double angle = atan2(col, row);
                    printf("Point (%d, %d): Phase = %.2f rad, %.2f°\n", row, col, angle, angle * (180 / M_PI));
                    symbol_index++;
                }
            }
        }
    }

    return 0;
}
