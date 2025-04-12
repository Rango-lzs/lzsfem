import os
import re
import argparse

def find_source_files(root_dir):
    """
    递归查找所有子目录中的 .h 和 .cpp 文件
    返回格式：保留子目录结构的相对路径 (使用正斜杠)
    """
    sources = []
    headers = []
    
    # 需要排除的目录 (如构建目录)
    exclude_dirs = {"build", "cmake-build", "third_party"}
    
    for dirpath, dirnames, filenames in os.walk(root_dir):
        # 过滤排除目录
        dirnames[:] = [d for d in dirnames if d not in exclude_dirs]
        
        for filename in filenames:
            # 生成相对于根目录的路径并统一为正斜杠
            full_path = os.path.join(dirpath, filename)
            rel_path = os.path.relpath(full_path, root_dir).replace('\\', '/')
            
            if filename.endswith(".cpp"):
                sources.append(rel_path)
            elif filename.endswith((".h", ".hpp")):
                headers.append(rel_path)
    
    return sorted(sources), sorted(headers)

def update_cmakelists(cmake_path, sources, headers):
    """
    更新 CMakeLists.txt 文件
    保留现有内容并智能更新变量
    """
    # 生成 CMake 格式的变量定义
    sources_str = "\n    ".join(f'"{s}"' for s in sources)
    headers_str = "\n    ".join(f'"{h}"' for h in headers)
    
    cmake_content = []
    if sources:
        cmake_content.append(f"set(SOURCES\n    {sources_str}\n)")
    if headers:
        cmake_content.append(f"set(HEADERS\n    {headers_str}\n)")
    
    # 读取原始文件内容
    try:
        with open(cmake_path, 'r', encoding='utf-8') as f:
            original = f.read()
    except FileNotFoundError:
        original = ""
    
    # 替换已有变量或追加新定义
    new_content = original
    if sources:
        new_content = re.sub(
            r'^set\(SOURCES\b.*?\)',
            cmake_content[0],
            new_content,
            flags=re.DOTALL | re.MULTILINE
        ) or new_content
        if "set(SOURCES" not in new_content:
            new_content = cmake_content[0] + "\n\n" + new_content
    
    if headers:
        new_content = re.sub(
            r'^set\(HEADERS\b.*?\)',
            cmake_content[1] if len(cmake_content)>1 else "",
            new_content,
            flags=re.DOTALL | re.MULTILINE
        ) or new_content
        if "set(HEADERS" not in new_content and headers:
            header_block = cmake_content[1] if len(cmake_content)>1 else ""
            new_content = header_block + "\n\n" + new_content
    
    # 确保 add_executable 包含所有文件
    new_content = re.sub(
        r'add_executable\(.*?\)',
        f'add_executable(${{PROJECT_NAME}}\n    ${{SOURCES}}\n    ${{HEADERS}}\n)',
        new_content,
        flags=re.DOTALL
    )
    
    # 写入更新后的内容
    with open(cmake_path, 'w', encoding='utf-8') as f:
        f.write(new_content.strip() + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='CMake源文件自动更新工具')
    parser.add_argument('--root', default='.', 
                       help='项目根目录 (默认当前目录)')
    parser.add_argument('--cmake', default='CMakeLists.txt',
                       help='CMakeLists.txt 路径 (默认当前目录)')
    
    args = parser.parse_args()
    
    # 规范化路径处理
    root_dir = os.path.abspath(args.root)
    cmake_path = os.path.join(root_dir, args.cmake)
    
    sources, headers = find_source_files(root_dir)
    print(f"扫描到 {len(sources)} 个源文件, {len(headers)} 个头文件")
    
    update_cmakelists(cmake_path, sources, headers)
    print(f"已更新 {cmake_path}")