# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
TensorFlow / Keras backend implementation.

All TensorFlow imports are deferred to method bodies.  The backend uses
the Keras functional API exclusively.

Input convention
----------------
pycsamt uses **channels-first** ordering internally (matching PyTorch).
This backend automatically transposes 1-D inputs

  ``(batch, channels, length)  →  (batch, length, channels)``

before passing them to Keras Conv1D layers, and transposes back if
needed.  The stored model handles this transparently.
"""
from __future__ import annotations

from typing import Any

import numpy as np

from ._base import BackendSpec, NeuralBackend
from ._detect import probe_backend

__all__ = ["TensorFlowBackend"]


class TensorFlowBackend(NeuralBackend):
    """TensorFlow/Keras backend for pycsamt AI/ML modules."""

    name = "tensorflow"

    @classmethod
    def is_available(cls) -> bool:
        return probe_backend("tensorflow")

    # ─── device ───────────────────────────────────────────────────────────

    def resolve_device(self, device: str | None = None) -> str:
        if device is not None:
            return device
        try:
            import tensorflow as tf
            gpus = tf.config.list_physical_devices("GPU")
            if gpus:
                return "/GPU:0"
        except ImportError:
            pass
        return "/CPU:0"

    # ─── model building ───────────────────────────────────────────────────

    def build(self, spec: BackendSpec) -> Any:
        self.require()
        arch = spec["arch"].lower()
        builders = {
            "cnn1d":    self._build_cnn1d,
            "cnn":      self._build_cnn1d,
            "resnet1d": self._build_resnet1d,
            "resnet":   self._build_resnet1d,
            "fcn1d":    self._build_fcn1d,
            "fcn":      self._build_fcn1d,
            "unet2d":   self._build_unet2d,
            "unet":     self._build_unet2d,
            "drcnn":    self._build_drcnn,
        }
        if arch not in builders:
            raise ValueError(
                f"Unknown arch {arch!r} for tensorflow backend.  "
                f"Available: {list(builders)}"
            )
        return builders[arch](spec)

    # ─── architecture implementations ─────────────────────────────────────

    def _build_cnn1d(self, spec) -> Any:
        """
        CNN1D — Puzyrev (2019/2021) style.

        Input shape: ``(n_features,)`` — flat feature vector.
        Internally reshaped to ``(n_features, 1)`` for Conv1D.
        """
        import tensorflow as tf
        from tensorflow.keras import Model, layers

        n_features = spec["n_features"]
        n_out = spec["n_out"]
        channels = spec.get("channels", (32, 64, 128))
        kernel_size = spec.get("kernel_size", 5)
        dropout = spec.get("dropout", 0.3)

        inp = tf.keras.Input(shape=(n_features,), name="input")
        x = layers.Reshape((n_features, 1))(inp)

        for ch in channels:
            x = layers.Conv1D(ch, kernel_size, padding="same")(x)
            x = layers.BatchNormalization()(x)
            x = layers.ReLU()(x)
            x = layers.MaxPooling1D(2, padding="same")(x)

        x = layers.Flatten()(x)
        x = layers.Dense(256, activation="relu")(x)
        x = layers.Dropout(dropout)(x)
        x = layers.Dense(128, activation="relu")(x)
        out = layers.Dense(n_out, name="output")(x)

        return Model(inp, out, name="cnn1d")

    def _build_resnet1d(self, spec) -> Any:
        """
        ResNet1D — Liu (2021 I8RFCN) style.

        Residual blocks with optional stride-2 downsampling and
        projection shortcuts.
        """
        import tensorflow as tf
        from tensorflow.keras import Model, layers

        n_features = spec["n_features"]
        n_out = spec["n_out"]
        channels = spec.get("channels", (64, 128, 256))
        dropout = spec.get("dropout", 0.3)
        n_blocks = spec.get("n_blocks", 2)

        def _res_block(x, out_ch, stride=1, name=""):
            shortcut = x
            x = layers.Conv1D(out_ch, 3, padding="same", strides=stride,
                              name=f"{name}_conv1")(x)
            x = layers.BatchNormalization(name=f"{name}_bn1")(x)
            x = layers.ReLU()(x)
            x = layers.Conv1D(out_ch, 3, padding="same",
                              name=f"{name}_conv2")(x)
            x = layers.BatchNormalization(name=f"{name}_bn2")(x)
            if shortcut.shape[-1] != out_ch or stride > 1:
                shortcut = layers.Conv1D(out_ch, 1, strides=stride,
                                         name=f"{name}_proj")(shortcut)
                shortcut = layers.BatchNormalization(
                    name=f"{name}_proj_bn")(shortcut)
            return layers.ReLU()(layers.Add()([x, shortcut]))

        inp = tf.keras.Input(shape=(n_features,), name="input")
        x = layers.Reshape((n_features, 1))(inp)

        # Stem
        x = layers.Conv1D(channels[0], 7, padding="same", name="stem_conv")(x)
        x = layers.BatchNormalization(name="stem_bn")(x)
        x = layers.ReLU()(x)

        for si, ch in enumerate(channels):
            for bi in range(n_blocks):
                stride = 2 if (bi == 0 and si > 0) else 1
                x = _res_block(x, ch, stride, name=f"s{si}b{bi}")

        x = layers.GlobalAveragePooling1D()(x)
        x = layers.Dropout(dropout)(x)
        out = layers.Dense(n_out, name="output")(x)

        return Model(inp, out, name="resnet1d")

    def _build_fcn1d(self, spec) -> Any:
        """
        FCN1D — Moghadas (2020) style — fully convolutional.
        """
        import tensorflow as tf
        from tensorflow.keras import Model, layers

        n_features = spec["n_features"]
        n_out = spec["n_out"]
        channels = spec.get("channels", (32, 64, 128, 64))
        dropout = spec.get("dropout", 0.2)

        inp = tf.keras.Input(shape=(n_features,), name="input")
        x = layers.Reshape((n_features, 1))(inp)

        for ch in channels:
            x = layers.Conv1D(ch, 3, padding="same")(x)
            x = layers.BatchNormalization()(x)
            x = layers.ReLU()(x)

        # Bottleneck 1×1 conv
        x = layers.Conv1D(channels[-1], 1, padding="same")(x)
        x = layers.GlobalAveragePooling1D()(x)
        x = layers.Dropout(dropout)(x)
        out = layers.Dense(n_out, name="output")(x)

        return Model(inp, out, name="fcn1d")

    def _build_unet2d(self, spec) -> Any:
        """
        2-D U-Net — Oh et al. (2019/2020) style.

        Input: ``(n_freqs, n_stations, n_in)`` — channels-last 2D.
        """
        import tensorflow as tf
        from tensorflow.keras import Model, layers

        n_in = spec.get("n_in", 4)
        n_out = spec.get("n_out", 1)
        channels = spec.get("channels", (32, 64, 128, 256, 512))
        spec.get("dropout", 0.2)

        def _conv2_block(x, ch, d):
            x = layers.Conv2D(ch, 3, padding="same", activation="relu")(x)
            x = layers.BatchNormalization()(x)
            x = layers.Conv2D(ch, 3, padding="same", activation="relu")(x)
            x = layers.BatchNormalization()(x)
            x = layers.SpatialDropout2D(d)(x)
            return x

        inp_shape = spec.get("input_shape", None)
        if inp_shape is None:
            inp_shape = (None, None, n_in)

        inp = tf.keras.Input(shape=inp_shape, name="input")
        ch = list(channels)
        n_stages = len(ch) - 1

        # Encoder
        skips = []
        x = inp
        for i in range(n_stages):
            x = _conv2_block(x, ch[i], 0.1)
            skips.append(x)
            x = layers.MaxPooling2D(2)(x)

        # Bridge
        x = _conv2_block(x, ch[-1], 0.2)

        # Decoder
        for i in range(n_stages - 1, -1, -1):
            x = layers.UpSampling2D(2)(x)
            x = layers.Concatenate()([x, skips[i]])
            x = _conv2_block(x, ch[i], 0.1)

        out = layers.Conv2D(n_out, 1, name="output")(x)
        return Model(inp, out, name="unet2d")

    def _build_drcnn(self, spec) -> Any:
        """
        DRCNN — Guo et al. (2021) style multi-modal fusion.

        Accepts a list of input tensors; one per modality.
        """
        import tensorflow as tf
        from tensorflow.keras import Model, layers

        n_features_list = spec["n_features_list"]
        n_out = spec["n_out"]
        growth_rate = spec.get("growth_rate", 32)
        n_layers = spec.get("n_layers", 6)
        hidden_dim = spec.get("hidden_dim", 256)
        dropout = spec.get("dropout", 0.2)

        def _dense_block(x, out_dim):
            accumulated = [x]
            for _ in range(n_layers):
                cat = layers.Concatenate()(accumulated) if len(accumulated) > 1 else accumulated[0]
                h = layers.Dense(growth_rate, activation="relu")(cat)
                h = layers.BatchNormalization()(h)
                h = layers.Dropout(dropout)(h)
                accumulated.append(h)
            cat_all = layers.Concatenate()(accumulated)
            out = layers.Dense(out_dim, activation="relu")(cat_all)
            out = layers.BatchNormalization()(out)
            # residual
            skip = layers.Dense(out_dim)(x)
            return layers.Add()([out, skip])

        inputs = [
            tf.keras.Input(shape=(nf,), name=f"input_{i}")
            for i, nf in enumerate(n_features_list)
        ]
        encoded = [_dense_block(inp, hidden_dim) for inp in inputs]

        fused = layers.Concatenate()(encoded) if len(encoded) > 1 else encoded[0]
        fused = _dense_block(fused, hidden_dim)

        out = layers.Dense(hidden_dim // 2, activation="relu")(fused)
        out = layers.Dropout(dropout)(out)
        out = layers.Dense(n_out, name="output")(out)

        return Model(inputs, out, name="drcnn")

    # ─── training ─────────────────────────────────────────────────────────

    def train(
        self,
        model,
        X_train, y_train,
        X_val, y_val,
        *,
        epochs=100, batch_size=256, lr=1e-3, patience=20,
        grad_clip=None, device=None, verbose=True, loss="mse",
    ):
        self.require()
        import tensorflow as tf

        dev = self.resolve_device(device)

        with tf.device(dev):
            loss_fn = self._get_loss_fn(loss)
            model.compile(
                optimizer=tf.keras.optimizers.Adam(learning_rate=lr),
                loss=loss_fn,
            )
            callbacks = [
                tf.keras.callbacks.EarlyStopping(
                    monitor="val_loss",
                    patience=patience,
                    restore_best_weights=True,
                    min_delta=1e-6,
                ),
                tf.keras.callbacks.ReduceLROnPlateau(
                    monitor="val_loss",
                    factor=0.5,
                    patience=max(5, patience // 3),
                    min_lr=1e-6,
                ),
            ]
            hist = model.fit(
                X_train, y_train,
                validation_data=(X_val, y_val),
                epochs=epochs,
                batch_size=batch_size,
                callbacks=callbacks,
                verbose=1 if verbose else 0,
            )

        return {
            "train_loss": hist.history["loss"],
            "val_loss": hist.history["val_loss"],
        }

    # ─── inference ────────────────────────────────────────────────────────

    def predict(self, model, X, *, device=None, batch_size=256):
        self.require()
        import tensorflow as tf
        dev = self.resolve_device(device)
        X = np.asarray(X, dtype=np.float32)
        with tf.device(dev):
            out = model.predict(X, batch_size=batch_size, verbose=0)
        return np.asarray(out)

    # ─── serialisation ────────────────────────────────────────────────────

    def get_weights(self, model) -> dict[str, np.ndarray]:
        # Store indexed weights to preserve order for set_weights
        weights = model.get_weights()
        return {f"layer_{i:04d}": w for i, w in enumerate(weights)}

    def set_weights(self, model, weights: dict[str, np.ndarray]) -> None:
        ordered = [weights[k] for k in sorted(weights)]
        model.set_weights(ordered)

    # ─── internal ─────────────────────────────────────────────────────────

    @staticmethod
    def _get_loss_fn(loss: str):
        if loss in ("mse", "masked_mse"):
            return "mse"
        if loss == "cross_entropy":
            return "sparse_categorical_crossentropy"
        return loss
