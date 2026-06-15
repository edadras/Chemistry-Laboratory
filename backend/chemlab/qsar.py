"""
qsar.py — Real QSAR / property models (machine learning).

Trains a scikit-learn model (RandomForest by default) on Morgan fingerprints
to predict a molecular property — either classification (e.g. BBB penetration,
active/inactive) or regression (e.g. pIC50). This is genuine ML, not a rule:
it learns structure→property from data, reports cross-validated metrics, and
persists for reuse.

Train from any CSV with columns `smiles` and `label`. A small bundled BBB
dataset lets a demo model train out-of-the-box; for real targets, plug in a
ChEMBL/PubChem export (see external_db.py) or your own CSV.
"""
from __future__ import annotations

import csv
import os
from dataclasses import dataclass, field

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from rdkit import DataStructs
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.model_selection import cross_val_predict, cross_val_score
from sklearn.metrics import (roc_auc_score, accuracy_score, r2_score,
                             mean_absolute_error)

from .molecule import Molecule

_DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
DEMO_DATASET = os.path.join(_DATA_DIR, "bbbp_demo.csv")
FP_BITS = 1024
FP_RADIUS = 2
_FPGEN = rdFingerprintGenerator.GetMorganGenerator(radius=FP_RADIUS, fpSize=FP_BITS)


def featurize(smiles: str) -> np.ndarray | None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    fp = _FPGEN.GetFingerprint(mol)
    arr = np.zeros((FP_BITS,), dtype=np.int8)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr


@dataclass
class QSARModel:
    name: str = "model"
    task: str = "classification"     # or "regression"
    model: object = None
    metrics: dict = field(default_factory=dict)
    n_train: int = 0

    # ----- training ----------------------------------------------------
    @classmethod
    def train_from_csv(cls, path: str, task: str = "classification",
                       name: str | None = None) -> "QSARModel":
        X, y = [], []
        with open(path, newline="", encoding="utf-8") as f:
            for row in csv.DictReader(f):
                feat = featurize(row["smiles"])
                if feat is None:
                    continue
                X.append(feat)
                y.append(float(row["label"]))
        if len(X) < 8:
            raise ValueError("داده‌ی کافی برای آموزش نیست (حداقل ۸ نمونه)")
        X = np.array(X)
        y = np.array(y)
        name = name or os.path.splitext(os.path.basename(path))[0]
        return cls._fit(X, y, task, name)

    @classmethod
    def _fit(cls, X, y, task, name) -> "QSARModel":
        if task == "classification":
            est = RandomForestClassifier(n_estimators=300, random_state=42,
                                         class_weight="balanced", n_jobs=-1)
            k = min(5, int(min(np.bincount(y.astype(int)))))
            k = max(2, k)
            try:
                proba = cross_val_predict(est, X, y, cv=k, method="predict_proba",
                                          n_jobs=-1)[:, 1]
                preds = (proba >= 0.5).astype(int)
                metrics = {"cv_auc": round(float(roc_auc_score(y, proba)), 3),
                           "cv_accuracy": round(float(accuracy_score(y, preds)), 3),
                           "cv_folds": k}
            except Exception as e:
                metrics = {"warning": f"اعتبارسنجی متقابل ناموفق: {e}"}
            est.fit(X, y)
        else:
            est = RandomForestRegressor(n_estimators=300, random_state=42, n_jobs=-1)
            k = 5
            try:
                preds = cross_val_predict(est, X, y, cv=k, n_jobs=-1)
                metrics = {"cv_r2": round(float(r2_score(y, preds)), 3),
                           "cv_mae": round(float(mean_absolute_error(y, preds)), 3),
                           "cv_folds": k}
            except Exception as e:
                metrics = {"warning": f"اعتبارسنجی متقابل ناموفق: {e}"}
            est.fit(X, y)
        return cls(name=name, task=task, model=est, metrics=metrics, n_train=len(X))

    # ----- inference ---------------------------------------------------
    def predict(self, smiles: str) -> dict:
        feat = featurize(smiles)
        if feat is None:
            return {"ok": False, "error_fa": "SMILES نامعتبر"}
        x = feat.reshape(1, -1)
        if self.task == "classification":
            proba = float(self.model.predict_proba(x)[0, 1])
            return {"ok": True, "task": "classification",
                    "probability_active": round(proba, 3),
                    "prediction": int(proba >= 0.5),
                    "score": round(proba, 3)}
        val = float(self.model.predict(x)[0])
        return {"ok": True, "task": "regression", "value": round(val, 3),
                "score": round(val, 3)}

    def to_dict(self) -> dict:
        return {"name": self.name, "task": self.task, "n_train": self.n_train,
                "metrics": self.metrics}


# Lazily-trained shared demo model (BBB penetration).
_DEMO_MODEL: QSARModel | None = None


def demo_model() -> QSARModel:
    global _DEMO_MODEL
    if _DEMO_MODEL is None:
        _DEMO_MODEL = QSARModel.train_from_csv(DEMO_DATASET, task="classification",
                                               name="BBB_penetration_demo")
    return _DEMO_MODEL
