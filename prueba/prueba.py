# app.py
import streamlit as st
import pandas as pd
from omics_analyzer import OmicsAnalyzer
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans
from io import StringIO

st.title("🔬 Análisis PCA y Clustering de expresión génica")

uploaded_file = st.file_uploader("📂 Sube un archivo CSV o TXT", type=["csv", "txt", "tsv"])

if uploaded_file:
    try:
        
        with open("temp_file.csv", "wb") as f:
            f.write(uploaded_file.getbuffer())
        st.success("✅ Datos cargados correctamente.")
        st.session_state.analizador = OmicsAnalyzer("temp_file.csv")

        analizador = st.session_state.analizador

        # Menú
        # Inicializar las claves de estado si no existen
        for key in ["show_pca", "show_clustering"]:
            if key not in st.session_state:
                st.session_state[key] = False

        col1, col2, col3 = st.columns([1, 1, 1])

        with col1:
            if st.button("PCA"):
                st.session_state.show_pca = True
                st.session_state.show_clustering = False

        with col2:
            if st.button("Análisis cluster"):
                st.session_state.show_clustering = True
                st.session_state.show_pca = False

        with col3:
            if st.button("Reset", type="primary"):
                for key in list(st.session_state.keys()):
                    if key != "analizador":
                        del st.session_state[key]
                
                st.session_state["show_pca"] = False
                st.session_state["show_clustering"] = False


        # PCA
        if st.session_state.show_pca: 
            fig, pca_df = analizador.run_pca()  
            
            st.subheader(":blue[Análisis de componentes principales]")
            st.write("Matriz de correlación")
            st.pyplot(analizador.correlation_matrix())

            st.write("Matriz de componentes principales")
            st.write(pca_df)
            st.write("Gráfico PC1 vs PC2")
            st.pyplot(fig)

            # varianza explicada y acumulada
            st.write("Variiaza explicada y acumulada por las CPs")
            st.pyplot(analizador.varianza_explicada())


        # Clustering
        if st.session_state.show_clustering:

            st.subheader(":blue[Análisis de Clustering]")

            option = st.selectbox(
                "Selecciona el tipo de análisis cluster:",
                (
                    "Selecciona una opción...",
                    "Clustering jerárquico - Método de Ward",
                    "Clustering no jerárquico - Método K-means"
                )
            )

            if option == "Clustering jerárquico - Método de Ward":
                if st.button("Ejecutar clustering jerárquico"):
                    fig = analizador.hierarchical_clustering()
                    st.pyplot(fig)

            if option == "Clustering no jerárquico - Método K-means":
                method = st.radio(
                    "Selecciona el método para determinar el número óptimo de clústeres:",
                    ["Método del codo", "Método de la silueta"]
                )
                method_code = 0 if method == "Método del codo" else 1

                
                modelo = KMeans(n_init=10, random_state=42)

                # Mostrar gráfico del método de selección de k
                fig = analizador.optimal_number_clusters(model=modelo, method=method_code)
                st.pyplot(fig)

                k = st.slider("Selecciona el número de clústeres para aplicar K-means", 2, 10, 3)
        
                if st.button("Ejecutar clustering K-means"):
                    st.subheader(f"🧩 Clustering K-means (k={k})")
                    conditions= ["LM","MT","PN","LM","MT","PN","LM","MT","PN","NM","NM","NM"]
                    fig = analizador.non_hierarchical_cluster(conditions=conditions,n_clusters=k)
                    st.pyplot(fig)


        


    except Exception as e:
        st.error(f"❌ Error al procesar el archivo: {e}")
