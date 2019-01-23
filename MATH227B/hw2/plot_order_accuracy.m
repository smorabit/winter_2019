function ps = plot_order_accuracy(f,a,b,n_list,boundary, iterations)

ps = zeros(1,length(n_list));
p_list = [];

j = 1;
for n = n_list
    p_ = zeros(1, iterations);
    for i = 1:iterations
        p_(i) = spline_order_accuracy(f,a,b,n,boundary,false);
    end
    ps(j) = mean(p_);
    p_list = [p_list; p_];
    j = j + 1;
    
end

plot(n_list(1)*ones(1,iterations), p_list(1,:), 'bd', 'MarkerFaceColor', 'b', 'MarkerSize', 3)
hold on
for i = 2:length(n_list)
    plot(n_list(i)*ones(1,iterations), p_list(i,:), 'bd', 'MarkerFaceColor', 'b', 'MarkerSize', 3)
end
plot(n_list, ps, 'LineWidth', 1.5)
xlabel("n")
ylabel("p")
grid on
hold off


size(p_list)
p_list(1,:)

end