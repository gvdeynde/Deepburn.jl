function brent(f, a::AbstractFloat, b::AbstractFloat; tol::AbstractFloat=sqrt(eps(a)), max_iter::Int=100)

    result = nothing
    fa = f(a)
    fb = f(b)
    d = zero(a)
    e = zero(a)

    if fa*fb > 0

        println("Root must be brackted")

    else

        c = b
        fc = fb

        iter = 0
        converged = false

        while iter < max_iter && !converged
            if (fb > 0 && fc > 0) || (fb < 0 && fc < 0)
                c = a
                fc = fa
                e = b-a
                d = e
            end

            if abs(fc) < abs(fb)
                a = b
                b = c
                c = a
                fa = fb
                fb = fc
                fc = fa
            end

            tolcvg = 2*eps(b)+tol/2
            xm = (c-b)/2

            if abs(xm) <= tolcvg || fb == zero(a)

                result = b
                converged = true

            else

                if abs(e) >= tolcvg && abs(fa) > abs(fb)
                    s= fb/fa

                    if (a==c)
                        p = 2*xm*s
                        q = 1-s
                    else
                        q = fa/fc
                        r = fb/fc
                        p = s*(2*xm*q*(q-r)-(b-a)*(r-1))
                        q = (q-1)*(r-1)*(s-1)
                    end

                    if p > 0
                        q = -q
                    end
                    p = abs(p)

                    m1 = 3*xm*q-abs(tolcvg*q)
                    m2 = abs(e*q)

                    if m1 < m2
                        m = m1
                    else
                        m = m2
                    end

                    if 2*p < m
                        e = d
                        d = p/q
                    else
                        d = xm
                        e = d
                    end
                else
                    d = xm
                    e = d
                end
                a = b
                fa = fb

                if abs(d) > tolcvg
                    b += d
                else
                    if xm > 0
                        b += tolcvg
                    else
                        b -= tolcvg
                    end
                end
                fb = f(b)
            end
        iter += 1
        end

        return result
    end
end
