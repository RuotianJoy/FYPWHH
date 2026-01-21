function [g, h] = get_coeffs(year)
if ~exist('pyversion', 'file')
    error('当前 MATLAB 未启用 Python 支持，无法调用 Python 版 igrf13。');
end

try
    mod = py.importlib.import_module('igrf13');
catch
    error(['未找到 Python 模块 "igrf13"。', newline, ...
        '请在 MATLAB 使用的 Python 环境中安装：pip install igrf13']);
end

py_result = mod.get_coeffs(year);
py_g = py_result{1};
py_h = py_result{2};

keys = cell(py.list(py_g.keys()));
n_values = zeros(1, numel(keys));
for i = 1:numel(keys)
    key = keys{i};
    n_values(i) = double(key{1});
end
n_max = max(n_values);

g = zeros(n_max + 1, n_max + 1);
h = zeros(n_max + 1, n_max + 1);

for i = 1:numel(keys)
    key = keys{i};
    n = double(key{1});
    m = double(key{2});
    g(n + 1, m + 1) = double(py_g{key});
    h(n + 1, m + 1) = double(py_h{key});
end
end

