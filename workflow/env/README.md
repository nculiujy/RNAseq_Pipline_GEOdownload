# workflow/env — 项目内置软件目录

将 Picard 的 JAR 包放置在此目录中，使项目自包含，避免依赖全局路径。

## 目录结构

```
workflow/env/
└── picard-2.18.2/
    └── picard.jar        ← 放置 picard.jar 文件
```

## 获取方式

```bash
# 从公共软件目录复制
cp /home/public_software_annotation/software/picard-2.18.2/picard.jar \
   workflow/env/picard-2.18.2/picard.jar

# 或从官方 GitHub 下载
mkdir -p workflow/env/picard-2.18.2
wget -O workflow/env/picard-2.18.2/picard.jar \
    https://github.com/broadinstitute/picard/releases/download/2.18.2/picard.jar
```

## 配置说明

`config/config.yaml` 中默认路径为：
```yaml
picard_jar: workflow/env/picard-2.18.2/picard.jar
```

若希望使用全局安装路径，修改为：
```yaml
picard_jar: /home/public_software_annotation/software/picard-2.18.2/picard.jar
```
