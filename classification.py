import scipy.io as sio
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
import time
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler

class PyTorchSVM:
    """使用PyTorch实现的SVM分类器"""

    def __init__(self, C=1.0, kernel='linear', max_iter=1000, lr=0.01, device='cuda'):
        """
        Args:
            C: 正则化参数
            kernel: 核函数类型 ('linear', 'rbf')
            max_iter: 最大迭代次数
            lr: 学习率
            device: 计算设备 ('cuda' 或 'cpu')
        """
        self.C = C
        self.kernel = kernel
        self.max_iter = max_iter
        self.lr = lr
        self.device = device if torch.cuda.is_available() else 'cpu'
        self.models = []  # 存储每个类别的二分类器
        self.scaler = StandardScaler()

    def rbf_kernel(self, X, Y, gamma=None):
        """RBF核函数"""
        if gamma is None:
            gamma = 1.0 / X.shape[1]  # 默认gamma值

        X_norm = torch.sum(X**2, dim=1).reshape(-1, 1)
        Y_norm = torch.sum(Y**2, dim=1).reshape(1, -1)
        K = torch.exp(-gamma * (X_norm + Y_norm - 2 * torch.mm(X, Y.t())))
        return K

    def fit(self, X, y):
        """训练SVM模型（使用One-vs-Rest策略）"""
        X = self.scaler.fit_transform(X)
        X_tensor = torch.FloatTensor(X).to(self.device)

        # 获取类别数量
        classes = np.unique(y)
        self.classes = classes

        print(f"训练 {len(classes)} 个SVM分类器 (One-vs-Rest)...")

        # 为每个类别训练一个二分类SVM
        for i, cls in enumerate(classes):
            print(f"训练类别 {cls+1}/{len(classes)}...")

            # 创建二分类标签
            y_binary = np.where(y == cls, 1, -1)
            y_tensor = torch.FloatTensor(y_binary).to(self.device)

            # 初始化权重和偏置
            n_samples, n_features = X_tensor.shape
            w = torch.zeros(n_features, device=self.device, requires_grad=True)
            b = torch.zeros(1, device=self.device, requires_grad=True)

            # 优化器
            optimizer = optim.Adam([w, b], lr=self.lr)

            # 训练循环
            for epoch in range(self.max_iter):
                optimizer.zero_grad()

                if self.kernel == 'linear':
                    # 线性核
                    outputs = torch.mm(X_tensor, w.view(-1, 1)).flatten() + b
                elif self.kernel == 'rbf':
                    # RBF核
                    K = self.rbf_kernel(X_tensor, X_tensor)
                    outputs = torch.mm(K, w.view(-1, 1)).flatten() + b

                # SVM合页损失
                loss = torch.mean(torch.clamp(1 - y_tensor * outputs, min=0)) + self.C * torch.norm(w)**2

                loss.backward()
                optimizer.step()

                if (epoch + 1) % 100 == 0:
                    print(f"  类别 {cls+1}, 迭代 {epoch+1}/{self.max_iter}, 损失: {loss.item():.4f}")

            # 保存模型参数
            self.models.append({
                'w': w.detach().clone(),
                'b': b.detach().clone(),
                'class': cls
            })

        print("所有SVM分类器训练完成!")
        return self

    def decision_function(self, X):
        """计算决策函数值"""
        X_scaled = self.scaler.transform(X)
        X_tensor = torch.FloatTensor(X_scaled).to(self.device)

        n_samples = X_tensor.shape[0]
        n_classes = len(self.models)
        decision_values = np.zeros((n_samples, n_classes))

        for i, model in enumerate(self.models):
            w, b, cls = model['w'], model['b'], model['class']

            if self.kernel == 'linear':
                values = torch.mm(X_tensor, w.view(-1, 1)).flatten() + b
            elif self.kernel == 'rbf':
                # 注意: 这里简化处理，实际应用中需要保存支持向量
                K = self.rbf_kernel(X_tensor, X_tensor)  # 这里需要修改以使用训练时的支持向量
                values = torch.mm(K, w.view(-1, 1)).flatten() + b

            decision_values[:, i] = values.cpu().numpy()

        return decision_values

    def predict(self, X):
        """预测类别"""
        decision_values = self.decision_function(X)
        predictions = np.argmax(decision_values, axis=1)
        return predictions

class PyTorchSpectralClassifier:
    """基于PyTorch SVM的高光谱图像分类器"""

    def __init__(self, spectral_lib_path, image_path, device='cuda'):
        self.device = device if torch.cuda.is_available() else 'cpu'
        print(f"使用设备: {self.device}")

        self.spectral_lib_path = spectral_lib_path
        self.image_path = image_path
        self.load_data()
        self.preprocess_data()

    def load_data(self):
        """加载MAT格式的光谱库和图像数据"""
        print("加载光谱库和图像数据...")

        # 加载光谱库
        spectral_data = sio.loadmat(self.spectral_lib_path)
        self.spectral_library = self._find_spectral_data(spectral_data)

        # 加载图像
        image_data = sio.loadmat(self.image_path)
        self.image_data = self._find_image_data(image_data)

        # 获取维度信息
        self.L, self.M = self.spectral_library.shape  # L:波段数, M:端元数
        self.L_img, self.N = self.image_data.shape    # L_img:图像波段数, N:像元数

        print(f"光谱库尺寸: {self.spectral_library.shape} (波段数 × 端元数)")
        print(f"图像数据尺寸: {self.image_data.shape} (波段数 × 像元数)")

        if self.L != self.L_img:
            raise ValueError("光谱库和图像的波段数不匹配！")

    def _find_spectral_data(self, data):
        """在mat文件中查找光谱库数据"""
        possible_keys = ['spectral_library', 'endmembers', 'library', 'spectra']
        for key in possible_keys:
            if key in data and isinstance(data[key], np.ndarray):
                return data[key]

        # 如果没找到标准名称，找第一个合适的二维数组
        for key in data:
            if (isinstance(data[key], np.ndarray) and
                    data[key].ndim == 2 and
                    not key.startswith('__')):
                print(f"使用变量 '{key}' 作为光谱库")
                return data[key]

        raise ValueError("在光谱库文件中找不到合适的光谱数据")

    def _find_image_data(self, data):
        """在mat文件中查找图像数据"""
        possible_keys = ['image', 'data', 'hyperspectral_image', 'img']
        for key in possible_keys:
            if key in data and isinstance(data[key], np.ndarray):
                return data[key]

        for key in data:
            if (isinstance(data[key], np.ndarray) and
                    data[key].ndim == 2 and
                    not key.startswith('__')):
                print(f"使用变量 '{key}' 作为图像数据")
                return data[key]

        raise ValueError("在图像文件中找不到合适的图像数据")

    def preprocess_data(self):
        """数据预处理"""
        print("数据预处理...")

        # 转置数据以适应机器学习格式 (样本数 × 特征数)
        self.X_train = self.spectral_library.T.astype(np.float32)  # M × L
        self.y_train = np.arange(self.M, dtype=np.int32)           # 0, 1, 2, ..., M-1

        # 图像数据转置
        self.X_image = self.image_data.T.astype(np.float32)        # N × L

        print(f"训练数据形状: {self.X_train.shape}")
        print(f"训练标签形状: {self.y_train.shape}")
        print(f"待分类图像数据形状: {self.X_image.shape}")

    def create_augmented_training_data(self, n_augmented=50, noise_level=0.01):
        """创建增强训练数据"""
        print(f"创建增强训练数据，每个类别 {n_augmented} 个样本...")

        augmented_features = []
        augmented_labels = []

        for class_idx in range(self.M):
            base_spectrum = self.X_train[class_idx]

            # 生成增强样本
            for i in range(n_augmented):
                # 添加高斯噪声
                noise = np.random.normal(0, noise_level, base_spectrum.shape)
                augmented_spectrum = base_spectrum + noise
                augmented_features.append(augmented_spectrum)
                augmented_labels.append(class_idx)

        self.X_train_augmented = np.array(augmented_features, dtype=np.float32)
        self.y_train_augmented = np.array(augmented_labels, dtype=np.int32)

        print(f"增强训练数据形状: {self.X_train_augmented.shape}")

        return self.X_train_augmented, self.y_train_augmented

    def train_svm(self, kernel='linear', C=1.0, use_augmented=True):
        """训练PyTorch SVM分类器"""
        print(f"训练PyTorch SVM分类器 (kernel={kernel}, C={C})...")

        if use_augmented and hasattr(self, 'X_train_augmented'):
            X_train = self.X_train_augmented
            y_train = self.y_train_augmented
        else:
            X_train = self.X_train
            y_train = self.y_train

        # 创建SVM分类器
        self.svm_classifier = PyTorchSVM(
            C=C,
            kernel=kernel,
            max_iter=500,
            lr=0.01,
            device=self.device
        )

        # 训练
        start_time = time.time()
        self.svm_classifier.fit(X_train, y_train)
        training_time = time.time() - start_time

        print(f"PyTorch SVM训练完成，耗时: {training_time:.2f}秒")

        self.training_time = training_time
        return self.svm_classifier

    def predict_image(self, batch_size=50000):
        """对整幅图像进行分类预测"""
        print("开始图像分类预测...")

        # 分批处理以避免内存问题
        n_pixels = self.X_image.shape[0]
        predictions = np.zeros(n_pixels, dtype=np.int32)

        start_time = time.time()

        for i in range(0, n_pixels, batch_size):
            end_idx = min(i + batch_size, n_pixels)
            batch = self.X_image[i:end_idx]

            # 预测当前批次
            batch_predictions = self.svm_classifier.predict(batch)
            predictions[i:end_idx] = batch_predictions

            # 进度显示
            progress = (end_idx / n_pixels) * 100
            if (i // batch_size) % 10 == 0:
                print(f"预测进度: {end_idx}/{n_pixels} ({progress:.1f}%)")

        prediction_time = time.time() - start_time

        print(f"图像分类完成，耗时: {prediction_time:.2f}秒")
        print(f"处理速度: {n_pixels/prediction_time:.0f} 像元/秒")

        self.predictions = predictions
        self.prediction_time = prediction_time

        return predictions

    def save_results(self, output_path, spatial_shape=None):
        """保存分类结果为MAT文件"""
        print(f"保存分类结果到 {output_path}...")

        # 确保预测结果是1-indexed（从1开始，按照光谱库顺序）
        classification_result = self.predictions + 1

        # 准备保存的数据
        save_data = {
            'classification_result': classification_result.astype(np.float32),
            'spectral_library': self.spectral_library,
            'class_labels': np.arange(1, self.M + 1),  # 1, 2, ..., M
            'classification_method': 'PyTorch_SVM',
            'svm_parameters': {
                'kernel': self.svm_classifier.kernel,
                'C': self.svm_classifier.C,
            },
            'performance_metrics': {
                'training_time': self.training_time,
                'prediction_time': self.prediction_time,
                'total_pixels': self.N
            },
            'image_dimensions': f'{self.L} x {self.N}'
        }

        # 如果提供了空间形状，重塑为2D图像
        if spatial_shape is not None:
            h, w = spatial_shape
            if h * w == len(classification_result):
                classification_map_2d = classification_result.reshape(h, w)
                save_data['classification_map_2d'] = classification_map_2d
                print(f"分类结果已重塑为2D: {classification_map_2d.shape}")
            else:
                print(f"警告: 空间形状 {spatial_shape} 与数据大小 {len(classification_result)} 不匹配")

        # 保存到MAT文件
        sio.savemat(output_path, save_data)
        print(f"结果已保存到: {output_path}")

        return save_data

    def evaluate_class_distribution(self):
        """评估类别分布"""
        if not hasattr(self, 'predictions'):
            print("请先进行预测")
            return

        unique, counts = np.unique(self.predictions, return_counts=True)

        print("\n类别分布统计:")
        print("=" * 40)
        for cls, count in zip(unique, counts):
            percentage = (count / len(self.predictions)) * 100
            print(f"类别 {cls+1}: {count:8d} 像元 ({percentage:6.2f}%)")

        return dict(zip(unique, counts))

    def visualize_results(self, spatial_shape=None, figsize=(15, 6)):
        """可视化分类结果"""
        if not hasattr(self, 'predictions'):
            print("请先进行预测")
            return

        if spatial_shape is None:
            print("未提供空间形状，无法可视化")
            return

        h, w = spatial_shape
        if h * w != len(self.predictions):
            print(f"空间形状 {spatial_shape} 与数据大小 {len(self.predictions)} 不匹配")
            return

        # 重塑为2D图像
        classification_map = (self.predictions + 1).reshape(h, w)

        # 创建可视化
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=figsize)

        # 1. 显示原始图像（RGB合成）
        if self.L >= 3:
            # 使用前三个波段创建RGB图像
            rgb_image = self.image_data[[0, 1, 2], :].T.reshape(h, w, 3)
            # 对比度拉伸
            p2, p98 = np.percentile(rgb_image, (2, 98))
            rgb_display = np.clip((rgb_image - p2) / (p98 - p2), 0, 1)
            ax1.imshow(rgb_display)
            ax1.set_title('原始图像 (RGB合成)')
        else:
            # 使用第一个波段显示灰度图像
            gray_image = self.image_data[0, :].reshape(h, w)
            ax1.imshow(gray_image, cmap='gray')
            ax1.set_title('原始图像 (第一波段)')
        ax1.axis('off')

        # 2. 显示分类结果
        im = ax2.imshow(classification_map, cmap='jet', vmin=1, vmax=self.M)
        ax2.set_title('PyTorch SVM分类结果')
        ax2.axis('off')
        plt.colorbar(im, ax=ax2, fraction=0.046, pad=0.04)

        # 3. 显示类别分布
        class_counts = np.bincount(self.predictions + 1, minlength=self.M + 2)[1:self.M + 1]
        classes = np.arange(1, self.M + 1)
        ax3.bar(classes, class_counts, color=plt.cm.jet(np.linspace(0, 1, self.M)))
        ax3.set_xlabel('类别索引')
        ax3.set_ylabel('像元数量')
        ax3.set_title('各类别像元数量分布')
        ax3.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig('pytorch_svm_classification_results.png', dpi=300, bbox_inches='tight')
        plt.show()

def main():
    """主函数"""
    # 文件路径
    spectral_lib_path = r'.\\spectral_library_urban.mat'
    image_path = r".\\Urban_R162.mat"
    output_path = r'.\\pytorch_svm_classification_result.mat'

    # 空间尺寸（根据您的图像实际尺寸修改）
    spatial_shape = (307, 307)  # 例如: 高度 × 宽度

    try:
        # 初始化分类器
        print("初始化PyTorch SVM分类器...")
        classifier = PyTorchSpectralClassifier(spectral_lib_path, image_path)

        # 创建增强训练数据（推荐）
        print("\n创建增强训练数据...")
        classifier.create_augmented_training_data(n_augmented=100, noise_level=0.02)

        # 训练SVM
        print("\n训练SVM分类器...")
        classifier.train_svm(
            kernel='linear',  # 可选择 'linear' 或 'rbf'
            C=1.0,           # 正则化参数
            use_augmented=True
        )

        # 图像分类
        print("\n开始图像分类...")
        predictions = classifier.predict_image(batch_size=50000)

        # 保存结果
        print("\n保存分类结果...")
        save_data = classifier.save_results(output_path, spatial_shape)

        # 评估类别分布
        print("\n分析类别分布...")
        class_distribution = classifier.evaluate_class_distribution()

        # 可视化结果
        print("\n生成可视化结果...")
        classifier.visualize_results(spatial_shape)

        # 性能总结
        print("\n" + "="*50)
        print("分类完成总结")
        print("="*50)
        print(f"总像元数: {classifier.N}")
        print(f"光谱库端元数: {classifier.M}")
        print(f"波段数: {classifier.L}")
        print(f"训练时间: {classifier.training_time:.2f}秒")
        print(f"预测时间: {classifier.prediction_time:.2f}秒")
        print(f"总处理时间: {classifier.training_time + classifier.prediction_time:.2f}秒")
        print(f"结果文件: {output_path}")

    except Exception as e:
        print(f"错误: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()