# Tratamiento de datos
# ==============================================================================
import numpy as np
import pandas as pd
#import statsmodels.api as sm
import os

# Gráficos
# ==============================================================================
import matplotlib.pyplot as plt
import matplotlib.font_manager
from matplotlib import style
style.use('ggplot') or plt.style.use('ggplot')
import seaborn as sns

# Preprocesado y modelado
# ==============================================================================
from sklearn.decomposition import PCA
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import scale
from sklearn.metrics import pairwise_distances
from scipy.cluster.hierarchy import linkage, dendrogram
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

# ==============================================================================

class OmicsAnalyzer:

    def __init__(self, filepath):
        self.data = self.load_data(filepath)
        self.normalized = None
        self.pca_result = None
        self.pca_scores = None

    def load_data(self, filepath):
        extension = os.path.splitext(filepath)[1]
        if extension == ".csv":
            df = pd.read_csv(filepath, sep=r'\s+', engine='python',index_col=0)
        elif extension in [".txt", ".tsv"]:
            df = pd.read_table(filepath, index_col=0)
        else:
            raise ValueError("Formato de archivo no soportado")
        return df.loc[df.sum(axis=1) > 0]  # quitar genes sin expresión
    
    def correlation_matrix(self):

        corr = self.data.corr()
        fig = sns.clustermap(corr, cmap="coolwarm", annot=False, method="ward")
        plt.title("Matriiz de correlación entre las muestras")

        return fig

    def normalizeLog2(self):
        cpm = self.data.div(self.data.sum(axis=0), axis=1) * 1e6
        log_cpm = np.log2(cpm + 1)
        self.normalized = log_cpm.T  # muestras como filas
        return self.normalized

    def run_pca(self, conditions=None):
        X = self.normalizeLog2()

        # PCA pipeline con escalado + reducción
        pca_pipe = make_pipeline(StandardScaler(), PCA())
        self.pca_result = pca_pipe.fit(X)
        self.pca_scores = self.pca_result.transform(X)


        # Preparar DataFrame con PCs
        n_components = self.pca_scores.shape[1]
        component_names = [f"PC{i+1}" for i in range(n_components)]
        pca_df = pd.DataFrame(self.pca_scores, columns=component_names, index=X.index)
        if conditions is not None and len(conditions) == len(X):
            pca_df["Grupo"] = conditions
            hue = "Grupo"
        else:
            pca_df["Muestra"] = X.index
            hue = "Muestra"

        # Graficar PCA
        fig, ax = plt.subplots(figsize=(8, 6))
        sns.set_theme(style="whitegrid")
        sns.scatterplot(data=pca_df, x="PC1", y="PC2", hue=hue, s=100, palette="Dark2", ax=ax)
        for i in range(len(pca_df)):
            ax.text(pca_df["PC1"][i]+0.5, pca_df["PC2"][i]+0.2, pca_df.index[i], fontsize=9)
        ax.set_title("PCA coloreado por muestra" if conditions is None else "PCA por grupo")
        ax.set_xlabel("PC1")
        ax.set_ylabel("PC2")
        plt.tight_layout()

        return fig, pca_df


    def varianza_explicada(self):

        # Varianza explicada
        pca = self.pca_result.named_steps['pca']
        var_exp = pca.explained_variance_ratio_
        cumulative_var_exp = np.cumsum(var_exp)

        # Crear figura
        fig, ax = plt.subplots(1, 2, figsize=(14, 5))

        # Scree plot
        ax[0].bar(range(1, len(var_exp) + 1), var_exp, color="skyblue")
        ax[0].set_title("Proporción de varianza explicada")
        ax[0].set_xlabel("Componentes principales")
        ax[0].set_ylabel("Proporción de varianza")
        ax[0].set_xticks(range(1, len(var_exp) + 1))
        ax[0].set_xticklabels([f"PC{i}" for i in range(1, len(var_exp)+1)], rotation=45)

        # Varianza acumulada
        ax[1].plot(range(1, len(cumulative_var_exp) + 1), cumulative_var_exp, marker='o', color="blue")
        ax[1].set_title("Varianza explicada acumulada")
        ax[1].set_xlabel("Número de Componentes")
        ax[1].set_ylabel("Varianza acumulada")
        ax[1].grid(True)

        plt.tight_layout()
        return fig


    def scale_data(self,axis=0):

        data = self.data

        if axis ==1:
            data = data.T
        
        datos_scaled = scale(X=data, axis=0, with_mean=True, with_std=True) 
        df_scaled = pd.DataFrame(datos_scaled, index=data.index, columns=data.columns)

        if axis == 1:
            df_scaled = df_scaled.T 

        return df_scaled
        
        


    def optimal_number_clusters (self,model,method):

        
        X = self.normalizeLog2()
        X = StandardScaler().fit_transform(X)

        k_range = range(2, 11)
        scores = []

        if method == 0:
            # Método del codo:
            for k in k_range:
                model.set_params(n_clusters=k)
                model.fit(X)
                scores.append(model.inertia_)

            fig, ax = plt.subplots()
            ax.plot(k_range, scores, marker='o')
            ax.set_xlabel("Número de clusters")
            ax.set_ylabel("Inercia (Suma de distancias intra-cluster)")
            ax.set_title("Método del codo")

        elif method == 1:
            # Método de la silueta
            for k in k_range:
                model.set_params(n_clusters=k)
                labels = model.fit_predict(X)
                score = silhouette_score(X, labels)
                scores.append(score)

            fig, ax = plt.subplots()
            ax.plot(k_range, scores, marker='o')
            ax.set_xlabel("Número de clusters")
            ax.set_ylabel("Score de silueta")
            ax.set_title("Método de la silueta")

        else:
            raise ValueError("El parámetro 'method' debe ser 0 (elbow) o 1 (silueta)")

        return fig


            
    def hierarchical_clustering(self):

        # Escalar
        X_scaled = self.normalizeLog2()
        X_scaled = StandardScaler().fit_transform(X_scaled)
    
        # Linkage para Ward
        Z = linkage(X_scaled, method='ward')

        # Dendrograma
        fig, ax = plt.subplots(figsize=(10, 5))
        dendrogram(Z, labels=self.data.T.index, leaf_rotation=90, ax=ax)
        ax.set_title("Clustering jerárquico (Ward)")
        ax.set_xlabel("Muestras")
        ax.set_ylabel("Distancia")
        return fig

    def non_hierarchical_cluster(self, conditions, n_clusters=3):
        # Escalamos las muestras (columnas) y transponemos para que cada muestra sea una fila
        
        X_scaled = self.normalizeLog2()
        X_scaled = StandardScaler().fit_transform(X_scaled)

        # K-means
        model = KMeans(n_clusters=n_clusters, n_init=10, random_state=42)
        labels = model.fit_predict(X_scaled)

        # Reducimos a 2D con PCA
        coords = PCA(n_components=2).fit_transform(X_scaled)

        # Usamos los nombres reales de muestra como etiquetas
        sample_names = list(self.data.columns)
        df = pd.DataFrame(coords, columns=["PC1", "PC2"], index=sample_names)
        df["Cluster"] = labels.astype(str)

        df["Tipo"] =  conditions

        # Gráfico
        fig, ax = plt.subplots()
        sns.scatterplot(data=df, x="PC1", y="PC2", style="Tipo", s=120, palette="Set2", ax=ax)

        # Añadir etiquetas con los nombres
        for i, row in df.iterrows():
            ax.text(row["PC1"] + 0.1, row["PC2"], i, fontsize=9)

        ax.set_title(f"K-means clustering (k={n_clusters})")
        return fig

    



