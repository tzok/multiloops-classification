import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split, cross_val_score, GridSearchCV
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, confusion_matrix, accuracy_score
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import roc_curve, auc
import warnings
warnings.filterwarnings('ignore')

class GraphClassificationPipeline:
    def __init__(self):
        self.models = {}
        self.scaler = StandardScaler()
        self.label_encoder = LabelEncoder()
        self.X_train = None
        self.X_test = None
        self.y_train = None
        self.y_test = None
        self.results = {}
        
    def load_data(self, file_path=None, X=None, y=None):
        """
        Wczytaj dane - albo z pliku albo z podanych macierzy
        """
        if file_path:
            data = pd.read_csv(file_path)
            # Zakładamy, że ostatnia kolumna to etykiety
            X = data.iloc[:, :-1]
            y = data.iloc[:, -1]
        
        print(f"Kształt danych: {X.shape}")
        print(f"Liczba klas: {len(np.unique(y))}")
        print(f"Rozkład klas:\n{pd.Series(y).value_counts()}")
        
        return X, y
    
    def preprocess_data(self, X, y, test_size=0.2, random_state=42):
        """
        Przygotowanie danych: podział train/test i normalizacja
        """
        # Kodowanie etykiet jeśli są stringami
        if y.dtype == 'object':
            y = self.label_encoder.fit_transform(y)
        
        # Podział na zbiory treningowy i testowy
        self.X_train, self.X_test, self.y_train, self.y_test = train_test_split(
            X, y, test_size=test_size, random_state=random_state, stratify=y
        )
        
        # Normalizacja cech
        self.X_train = self.scaler.fit_transform(self.X_train)
        self.X_test = self.scaler.transform(self.X_test)
        
        print(f"Rozmiar zbioru treningowego: {self.X_train.shape}")
        print(f"Rozmiar zbioru testowego: {self.X_test.shape}")
    
    def initialize_models(self):
        """
        Inicjalizacja klasyfikatorów zgodnie z diagramem
        """
        self.models = {
            'k-NN': KNeighborsClassifier(),
            'SVM': SVC(probability=True),
            'Random Forest': RandomForestClassifier(random_state=42)
        }
    
    def hyperparameter_tuning(self):
        """
        Dostrajanie hiperparametrów dla każdego modelu
        """
        param_grids = {
            'k-NN': {
                'n_neighbors': [3, 5, 7, 9, 11],
                'weights': ['uniform', 'distance'],
                'metric': ['euclidean', 'manhattan']
            },
            'SVM': {
                'C': [0.1, 1, 10, 100],
                'kernel': ['rbf', 'linear', 'poly'],
                'gamma': ['scale', 'auto']
            },
            'Random Forest': {
                'n_estimators': [50, 100, 200],
                'max_depth': [None, 10, 20, 30],
                'min_samples_split': [2, 5, 10]
            }
        }
        
        best_models = {}
        
        for name, model in self.models.items():
            print(f"Dostrajanie hiperparametrów dla {name}...")
            
            grid_search = GridSearchCV(
                model, param_grids[name], 
                cv=5, scoring='accuracy', 
                n_jobs=-1, verbose=0
            )
            
            grid_search.fit(self.X_train, self.y_train)
            best_models[name] = grid_search.best_estimator_
            
            print(f"Najlepsze parametry dla {name}: {grid_search.best_params_}")
            print(f"Najlepszy wynik CV: {grid_search.best_score_:.4f}")
            print("-" * 50)
        
        self.models = best_models
    
    def train_models(self):
        """
        Trenowanie wszystkich modeli
        """
        for name, model in self.models.items():
            print(f"Trenowanie {name}...")
            model.fit(self.X_train, self.y_train)
    
    def evaluate_models(self):
        """
        Ewaluacja modeli i zapisanie wyników
        """
        for name, model in self.models.items():
            # Predykcje
            y_pred = model.predict(self.X_test)
            y_prob = model.predict_proba(self.X_test) if hasattr(model, 'predict_proba') else None
            
            # Metryki
            accuracy = accuracy_score(self.y_test, y_pred)
            
            # Cross-validation score
            cv_scores = cross_val_score(model, self.X_train, self.y_train, cv=5)
            
            self.results[name] = {
                'accuracy': accuracy,
                'cv_mean': cv_scores.mean(),
                'cv_std': cv_scores.std(),
                'predictions': y_pred,
                'probabilities': y_prob,
                'classification_report': classification_report(self.y_test, y_pred)
            }
            
            print(f"\n{name} - Wyniki:")
            print(f"Dokładność na zbiorze testowym: {accuracy:.4f}")
            print(f"Cross-validation: {cv_scores.mean():.4f} (+/- {cv_scores.std() * 2:.4f})")
            print(f"Raport klasyfikacji:\n{self.results[name]['classification_report']}")
            print("-" * 80)
    
    def plot_results(self):
        """
        Wizualizacja wyników
        """
        # 1. Porównanie dokładności modeli
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        
        # Dokładność modeli
        model_names = list(self.results.keys())
        accuracies = [self.results[name]['accuracy'] for name in model_names]
        cv_means = [self.results[name]['cv_mean'] for name in model_names]
        
        axes[0, 0].bar(model_names, accuracies, alpha=0.7, label='Test Accuracy')
        axes[0, 0].bar(model_names, cv_means, alpha=0.7, label='CV Mean')
        axes[0, 0].set_title('Porównanie Dokładności Modeli')
        axes[0, 0].set_ylabel('Dokładność')
        axes[0, 0].legend()
        axes[0, 0].tick_params(axis='x', rotation=45)
        
        # Macierz konfuzji dla najlepszego modelu
        best_model_name = max(self.results.keys(), key=lambda x: self.results[x]['accuracy'])
        best_predictions = self.results[best_model_name]['predictions']
        
        cm = confusion_matrix(self.y_test, best_predictions)
        sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', ax=axes[0, 1])
        axes[0, 1].set_title(f'Macierz Konfuzji - {best_model_name}')
        axes[0, 1].set_ylabel('Prawdziwe etykiety')
        axes[0, 1].set_xlabel('Predykcje')
        
        # Cross-validation scores
        cv_data = []
        for name in model_names:
            model = self.models[name]
            cv_scores = cross_val_score(model, self.X_train, self.y_train, cv=5)
            cv_data.extend([(name, score) for score in cv_scores])
        
        cv_df = pd.DataFrame(cv_data, columns=['Model', 'CV Score'])
        sns.boxplot(data=cv_df, x='Model', y='CV Score', ax=axes[1, 0])
        axes[1, 0].set_title('Rozkład Wyników Cross-Validation')
        axes[1, 0].tick_params(axis='x', rotation=45)
        
        # ROC curves (dla klasyfikacji binarnej lub multiclass)
        n_classes = len(np.unique(self.y_test))
        if n_classes == 2:
            # Klasyfikacja binarna
            for name, model in self.models.items():
                if hasattr(model, 'predict_proba'):
                    y_score = model.predict_proba(self.X_test)[:, 1]
                    fpr, tpr, _ = roc_curve(self.y_test, y_score)
                    roc_auc = auc(fpr, tpr)
                    axes[1, 1].plot(fpr, tpr, label=f'{name} (AUC = {roc_auc:.2f})')
            
            axes[1, 1].plot([0, 1], [0, 1], 'k--')
            axes[1, 1].set_title('Krzywe ROC')
            axes[1, 1].set_xlabel('False Positive Rate')
            axes[1, 1].set_ylabel('True Positive Rate')
            axes[1, 1].legend()
        else:
            # Dla klasyfikacji wieloklasowej - pokazujemy feature importance dla Random Forest
            if 'Random Forest' in self.models:
                rf_model = self.models['Random Forest']
                feature_importance = rf_model.feature_importances_
                feature_names = [f'Feature {i+1}' for i in range(len(feature_importance))]
                
                axes[1, 1].bar(feature_names, feature_importance)
                axes[1, 1].set_title('Ważność Cech (Random Forest)')
                axes[1, 1].set_ylabel('Ważność')
                axes[1, 1].tick_params(axis='x', rotation=45)
        
        plt.tight_layout()
        plt.show()
    
    def get_best_model(self):
        """
        Zwraca najlepszy model na podstawie dokładności
        """
        best_model_name = max(self.results.keys(), key=lambda x: self.results[x]['accuracy'])
        return best_model_name, self.models[best_model_name]
    
    def run_pipeline(self, X, y, tune_hyperparameters=True):
        """
        Uruchomienie całego pipeline'u
        """
        print("=" * 80)
        print("PIPELINE KLASYFIKACJI GRAFÓW")
        print("=" * 80)
        
        # Przygotowanie danych
        print("\n1. PRZYGOTOWANIE DANYCH")
        print("-" * 40)
        self.preprocess_data(X, y)
        
        # Inicjalizacja modeli
        print("\n2. INICJALIZACJA MODELI")
        print("-" * 40)
        self.initialize_models()
        
        # Dostrajanie hiperparametrów (opcjonalne)
        if tune_hyperparameters:
            print("\n3. DOSTRAJANIE HIPERPARAMETRÓW")
            print("-" * 40)
            self.hyperparameter_tuning()
        
        # Trenowanie
        print("\n4. TRENOWANIE MODELI")
        print("-" * 40)
        self.train_models()
        
        # Ewaluacja
        print("\n5. EWALUACJA MODELI")
        print("-" * 40)
        self.evaluate_models()
        
        # Wizualizacja
        print("\n6. WIZUALIZACJA WYNIKÓW")
        print("-" * 40)
        self.plot_results()
        
        # Najlepszy model
        best_name, best_model = self.get_best_model()
        print(f"\n7. NAJLEPSZY MODEL: {best_name}")
        print(f"Dokładność: {self.results[best_name]['accuracy']:.4f}")
        
        return best_model

# Przykład użycia
if __name__ == "__main__":
    dataset = pd.read_csv('/home/matrament/projects/master_thesis/data/3way_junctions_features_graphs.csv')  # Wczytaj swój zbiór danych
    X = dataset.drop(columns=['family'])  
    y = dataset['family'] 
    
    # Uruchomienie pipeline'u
    pipeline = GraphClassificationPipeline()
    X_data, y_data = pipeline.load_data(X=X, y=y)
    best_model = pipeline.run_pipeline(X_data, y_data, tune_hyperparameters=True)
    
    # Predykcja na nowych danych
    # new_graph_features = np.array([[1.2, -0.5, 2.1, 0.8, -1.0, 1.5]])
    # prediction = best_model.predict(pipeline.scaler.transform(new_graph_features))
    # print(f"Predykcja dla nowego grafu: {prediction[0]}")