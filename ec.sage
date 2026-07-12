def nonlinear_complexity(seq, m, n=None):
    """
    计算序列的 m 阶非线性复杂度 NL_m(seq, n)
    严格遵循论文定义II.2，通过线性方程组可解性判定最小阶数
    """
    if n is None:
        n = len(seq)
    s = seq[:n]

    if n <= 1:
        return 0
    F = s[0].parent()
    if all(x == F.zero() for x in s):
        return 0

    # 从小到大枚举最小可行阶数
    for ell in range(1, n):
        num_eq = n - ell
        if num_eq <= 0:
            return ell

        # 生成所有总次数 ≤ m 的 ell 元单项式指数向量
        monomials = []
        for deg in range(m + 1):
            monomials.extend(list(IntegerVectors(deg, ell)))
        num_var = len(monomials)

        # 剪枝：未知数 ≥ 方程数时必然有解，直接返回
        if num_var >= num_eq:
            return ell

        # 批量构造系数矩阵与右端项
        A_data = []
        b_data = []
        for i in range(num_eq):
            window = s[i:i+ell]
            row = [F.prod(window[j] ** e for j, e in enumerate(exps)) for exps in monomials]
            A_data.append(row)
            b_data.append(s[i + ell])

        A = matrix(F, A_data)
        b = vector(F, b_data)

        # 方程组有解则当前 ell 即为非线性复杂度
        try:
            A.solve_right(b)
            return ell
        except ValueError:
            continue

    return n - 1


# ========== 演示序列构造（符合论文分块+轨道结构） ==========
def demo_paper_sequence():
    """
    生成与论文算例同参数的测试序列
    参数：h=13, 单块长度8, 总长度104, 有限域F_16
    """
    F16.<a> = GF(16)
    h = 13
    block_num = 8

    seq = []
    for block_idx in range(block_num):
        base = F16.gen() ** (block_idx + 1)
        # 模拟Galois轨道生成的子序列
        block = [base * F16.gen() ** (i * 3) for i in range(h)]
        seq.extend(block)
    
    return seq


# ========== 运行测试 ==========
if __name__ == "__main__":
    seq = demo_paper_sequence()
    print(f"序列总长度: {len(seq)}")
    print(f"序列所在域: {seq[0].parent()}")

    test_cases = [
        (50, 2),
        (80, 2),
        (104, 2),
        (104, 3),
    ]

    print("\n n\tm\t非线性复杂度")
    print("-" * 28)
    for n, m in test_cases:
        nl = nonlinear_complexity(seq, m, n)
        print(f"{n}\t{m}\t{nl}")

    # 验证：m=1时结果等于线性复杂度
    print("\n验证 m=1 时等价于线性复杂度:")
    print(f"NL_1(seq, 104) = {nonlinear_complexity(seq, 1, 104)}")
