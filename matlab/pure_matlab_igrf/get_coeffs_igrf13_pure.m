function [g, h] = get_coeffs_igrf13_pure(year)
persistent coeffs
if isempty(coeffs)
    coeffs = load_coeffs();
end
y = max(min(year, 2025.0), 1900.0);
years = coeffs.years;
K = numel(years);
n_max = coeffs.n_max;
gmat = zeros(n_max + 1, n_max + 1);
hmat = zeros(n_max + 1, n_max + 1);
if y <= years(1)
    gmat = coeffs.g(:, :, 1);
    hmat = coeffs.h(:, :, 1);
elseif y <= years(end)
    idx = find(years <= y, 1, 'last');
    if idx == K
        gmat = coeffs.g(:, :, K);
        hmat = coeffs.h(:, :, K);
    else
        t0 = years(idx);
        t1 = years(idx + 1);
        f = (y - t0) / (t1 - t0);
        g0 = coeffs.g(:, :, idx);
        g1 = coeffs.g(:, :, idx + 1);
        h0 = coeffs.h(:, :, idx);
        h1 = coeffs.h(:, :, idx + 1);
        gmat = g0 + f * (g1 - g0);
        hmat = h0 + f * (h1 - h0);
    end
else
    dt = y - years(end);
    if dt > 5
        dt = 5;
    end
    gmat = coeffs.g(:, :, K) + dt * coeffs.g_sv;
    hmat = coeffs.h(:, :, K) + dt * coeffs.h_sv;
end
g = gmat;
h = hmat;
end

function coeffs = load_coeffs()
fname = fullfile(fileparts(mfilename('fullpath')), 'igrf13coeffs.txt');
fid = fopen(fname, 'r');
if fid < 0
    error('未找到igrf13coeffs.txt，请将官方IGRF13系数文件放在pure_matlab_igrf文件夹下。');
end
lines = {};
while true
    line = fgetl(fid);
    if ~ischar(line)
        break
    end
    lines{end + 1} = line;
end
fclose(fid);
header_idx = 0;
for i = 1:numel(lines)
    line = strtrim(lines{i});
    if numel(line) >= 3 && strcmp(line(1:3), 'g/h')
        header_idx = i;
        break
    end
end
if header_idx == 0
    error('igrf13coeffs.txt格式不正确，未找到g/h表头行。');
end
header_tokens = strsplit(strtrim(lines{header_idx}));
labels = header_tokens(4:end);
num_labels = numel(labels);
years = zeros(1, num_labels - 1);
for j = 1:num_labels - 1
    years(j) = str2double(labels{j});
end
n_max = 0;
for i = header_idx + 1:numel(lines)
    line = strtrim(lines{i});
    if isempty(line)
        continue
    end
    if line(1) ~= 'g' && line(1) ~= 'h'
        continue
    end
    tokens = strsplit(line);
    if numel(tokens) < 4
        continue
    end
    n = str2double(tokens{2});
    if ~isnan(n) && n > n_max
        n_max = n;
    end
end
K = numel(years);
g = zeros(n_max + 1, n_max + 1, K);
h = zeros(n_max + 1, n_max + 1, K);
g_sv = zeros(n_max + 1, n_max + 1);
h_sv = zeros(n_max + 1, n_max + 1);
for i = header_idx + 1:numel(lines)
    line = strtrim(lines{i});
    if isempty(line)
        continue
    end
    if line(1) ~= 'g' && line(1) ~= 'h'
        continue
    end
    tokens = strsplit(line);
    if numel(tokens) < 4 + K
        continue
    end
    ttype = tokens{1};
    n = str2double(tokens{2});
    m = str2double(tokens{3});
    ni = n + 1;
    mi = m + 1;
    for j = 1:K
        vstr = tokens[3 + j];
        if strcmp(vstr, '-')
            v = 0;
        else
            v = str2double(vstr);
        end
        if ttype == 'g'
            g(ni, mi, j) = v;
        else
            h(ni, mi, j) = v;
        end
    end
    sv_str = tokens{end};
    if strcmp(sv_str, '-')
        sv_val = 0;
    else
        sv_val = str2double(sv_str);
    end
    if ttype == 'g'
        g_sv(ni, mi) = sv_val;
    else
        h_sv(ni, mi) = sv_val;
    end
end
coeffs.years = years;
coeffs.n_max = n_max;
coeffs.g = g;
coeffs.h = h;
coeffs.g_sv = g_sv;
coeffs.h_sv = h_sv;
end

