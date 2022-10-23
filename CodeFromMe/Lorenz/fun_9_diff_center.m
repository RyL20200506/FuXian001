% 做中心差分的
function diff_center = fun_9_diff_center(known_data)
diff_before = known_data(2:end-1) - known_data(1:end-2);  % 前向差分: #验证无误
diff_after = (known_data(3:end) - known_data(2:end-1));  % 后向差分  #验证无误
diff_center = (diff_before+diff_after)/2;  % 中心差分
end
