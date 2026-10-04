def load_envi_data(self):
    # ハイパースペクトルデータを読み込みと前処理（正規化・輝度リセット）
    try:
        # ENVIファイルの読み込み、float32へのキャスト、形状検証、不正値置換
        img = envi.open(file_path)
        raw_data = np.array(img.load(), dtype=np.float32)

        # 0.0 ~ 1.0 の範囲に正規化
        d_min, d_max = np.min(raw_data), np.max(raw_data)
        if d_max > d_min:
            self.hsi_data = (raw_data - d_min) / (d_max - d_min)
        else:
            self.hsi_data = raw_data

        self.height, self.width, self.bands = self.hsi_data.shape

        # 初期擬似RGBのバンド設定
        self.r_idx = min(self.bands - 1, max(0, int(self.bands * 0.8)))
        self.g_idx = min(self.bands - 1, max(0, int(self.bands * 0.5)))
        self.b_idx = min(self.bands - 1, max(0, int(self.bands * 0.2)))

def on_select_roi(self, eclick, erelease):
    # ROI選択時に呼び出されるコールバック
    # 選択範囲の(x1, y1)・(x2, y2)を正規化し、画像サイズ内に制限
    x1, y1 = int(eclick.xdata), int(eclick.ydata)
    x2, y2 = int(erelease.xdata), int(erelease.ydata)
    x1, x2 = min(x1, x2), max(x1, x2)
    y1, y2 = min(y1, y2), max(y1, y2)
    x1, x2 = max(0, x1), min(self.width, x2)
    y1, y2 = max(0, y1), min(self.height, y2)

    # 選択領域内の平均スペクトルを抽出・算出・表示
    roi_pixels = self.hsi_data[y1:y2, x1:x2, :]
    raw_specs = roi_pixels.reshape(-1, self.bands)
    mean_spec = np.mean(raw_specs, axis=0)

def run_segmentation(self, thresh):
    # セグメンテーション実行部
    try:
        method_idx = self.combo_method.current()
        h, w, b = self.height, self.width, self.bands
        flat_data = self.hsi_data.reshape(-1, b)

        result_map = np.zeros(h * w, dtype=np.int32)

        # 大容量データでのメモリ不足（OOM）防止と進捗更新のためのチャンク分割サイズ設定
        chunk_size = max(1000, (h * w) // 100)

        mean_specs = np.array([roi['mean_spec'] for roi in self.rois])

        if method_idx == 0:  # ユークリッド距離
            for i in range(0, flat_data.shape[0], chunk_size):
                chunk = flat_data[i:i + chunk_size]
                dists = np.linalg.norm(chunk[:, np.newaxis, :] - mean_specs[np.newaxis, :, :], axis=2)
                min_dists = np.min(dists, axis=1)
                best_roi = np.argmin(dists, axis=1) + 1

                mask = min_dists <= thresh
                chunk_res = np.zeros(len(chunk), dtype=np.int32)
                chunk_res[mask] = best_roi[mask]
                result_map[i:i + chunk_size] = chunk_res
                update_progress(min(100, ((i + len(chunk)) / (h * w)) * 100))

        elif method_idx == 1:  # スペクトル角（SAM）
            norm_means = np.linalg.norm(mean_specs, axis=1)
            norm_means[norm_means == 0] = 1e-10

            for i in range(0, flat_data.shape[0], chunk_size):
                chunk = flat_data[i:i + chunk_size]
                norm_chunk = np.linalg.norm(chunk, axis=1)
                norm_chunk[norm_chunk == 0] = 1e-10

                dot_prod = np.dot(chunk, mean_specs.T)
                denom = norm_chunk[:, np.newaxis] * norm_means[np.newaxis, :]
                cos_angle = np.clip(dot_prod / denom, -1.0, 1.0)
                angles = np.arccos(cos_angle)

                min_angles = np.min(angles, axis=1)
                best_roi = np.argmin(angles, axis=1) + 1

                mask = min_angles <= thresh
                chunk_res = np.zeros(len(chunk), dtype=np.int32)
                chunk_res[mask] = best_roi[mask]
                result_map[i:i + chunk_size] = chunk_res
                update_progress(min(100, ((i + len(chunk)) / (h * w)) * 100))

        elif method_idx == 2:  # マハラノビス距離
            roi_params = []
            for roi in self.rois:
                raw_specs = roi['raw_specs']
                mean_spec = roi['mean_spec']

                # 各ROIの共分散行列計算
                # 特異行列による逆行列計算エラー回避のため、sqrt(0.01) = 0.1 (10%相当) の許容ノイズマージンを各バンドに付与
                # 反射率スケール(0.0〜1.0)において直観的なスケール感と一致させている
                cov = np.cov(raw_specs, rowvar=False)
                cov += np.eye(self.bands) * 1e-2
                inv_cov = np.linalg.inv(cov)

                roi_params.append((mean_spec, inv_cov))

            total_pixels = flat_data.shape[0]

            for i in range(0, total_pixels, chunk_size):
                chunk = flat_data[i:i + chunk_size]
                chunk_min_dists = np.full(chunk.shape[0], np.inf)
                chunk_result = np.zeros(chunk.shape[0], dtype=int)

                for idx, (mean_spec, inv_cov) in enumerate(roi_params):
                    diff = chunk - mean_spec
                    dist_sq = np.sum(np.dot(diff, inv_cov) * diff, axis=1)
                    dist = np.sqrt(np.maximum(0, dist_sq) / self.bands)

                    mask = (dist < chunk_min_dists) & (dist < thresh)
                    chunk_min_dists[mask] = dist[mask]
                    chunk_result[mask] = idx + 1

                result_map[i:i + chunk_size] = chunk_result
                update_progress(min(100, ((i + len(chunk)) / total_pixels) * 100))

        elif method_idx == 3:  # 線形SVM
            X_train, y_train = [], []
            # 計算速度維持のため、学習データをROIごとに上限設定してランダムサンプリング
            max_samples_per_roi = 10000

            for idx, roi in enumerate(self.rois):
                raw_specs = roi['raw_specs']
                n_samples = len(raw_specs)

                if n_samples > max_samples_per_roi:
                    indices = np.random.choice(n_samples, max_samples_per_roi, replace=False)
                    sampled_specs = raw_specs[indices]
                else:
                    sampled_specs = raw_specs

                X_train.append(sampled_specs)
                y_train.append(np.full(len(sampled_specs), idx + 1))

            X_train = np.vstack(X_train)
            y_train = np.concatenate(y_train)

            clf = LinearSVC(dual=False, max_iter=2000)
            clf.fit(X_train, y_train)

            for i in range(0, flat_data.shape[0], chunk_size):
                chunk = flat_data[i:i + chunk_size]
                result_map[i:i + chunk_size] = clf.predict(chunk)
                update_progress(min(100, ((i + len(chunk)) / (h * w)) * 100))

        elif method_idx == 4:  # 非線形SVM
            X_train, y_train = [], []
            # 計算速度維持のため、学習データをROIごとに上限設定してランダムサンプリング
            max_samples_per_roi = 1000

            for idx, roi in enumerate(self.rois):
                raw_specs = roi['raw_specs']
                n_samples = len(raw_specs)

                if n_samples > max_samples_per_roi:
                    indices = np.random.choice(n_samples, max_samples_per_roi, replace=False)
                    sampled_specs = raw_specs[indices]
                else:
                    sampled_specs = raw_specs

                X_train.append(sampled_specs)
                y_train.append(np.full(len(sampled_specs), idx + 1))

            X_train = np.vstack(X_train)
            y_train = np.concatenate(y_train)

            clf = SVC(kernel='rbf', gamma='scale')
            clf.fit(X_train, y_train)

            for i in range(0, flat_data.shape[0], chunk_size):
                chunk = flat_data[i:i + chunk_size]
                result_map[i:i + chunk_size] = clf.predict(chunk)
                update_progress(min(100, ((i + len(chunk)) / (h * w)) * 100))

        update_progress(100)
        seg_result = result_map.reshape((h, w))

        seg_rgb = np.zeros((h, w, 3), dtype=np.float32)
        for idx, roi in enumerate(self.rois):
            rgb_color = mcolors.to_rgb(roi['color'])
            seg_rgb[seg_result == (idx + 1)] = rgb_color

        self.root.after(0, self.show_segmentation_result, seg_rgb)

def show_segmentation_result(self, seg_rgb):
    # 結果画像とボタンの描画設定
    fig_res, ax_res = plt.subplots(figsize=(8, 6))
    fig_res.subplots_adjust(left=0, right=1, top=1, bottom=0)
    ax_res.imshow(seg_rgb)
    ax_res.axis("off")

    btn_frame = tk.Frame(res_win, width=btn_w, height=btn_h)
    btn_frame.pack_propagate(False)
    btn_frame.pack(side=tk.BOTTOM, pady=5)

    def save_result():
        # 結果画像のファイル出力
        save_path = filedialog.asksaveasfilename(
            parent=res_win,
            defaultextension=".png", 
            filetypes=[("PNG file", "*.png"), ("JPEG file", "*.jpg"), ("BMP file", "*.bmp")]
        )
        if save_path:
            try:
                img_uint8 = (seg_rgb * 255).astype(np.uint8)
                img_pil = Image.fromarray(img_uint8)
                img_pil.save(save_path)
                messagebox.showinfo("Saved", f"結果画像を保存しました:\n{save_path}", parent=res_win)
            except Exception as e:
                messagebox.showerror("Error", f"画像の保存に失敗しました:\n{e}", parent=res_win)
