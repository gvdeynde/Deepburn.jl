function brentzero(f, a::AbstractFloat, b::AbstractFloat; tol::AbstractFloat=sqrt(eps(a)), max_iter::Int=100)

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

function brentmin(f, a::AbstractFloat, b::AbstractFloat; tol::AbstractFloat=sqrt(eps(a)), max_iter::Int=100)

    if a > b
        a, b = b, a
    end

    c = (3 - sqrt(5))/2
    d = 0

    sa = a
    sb = b
    x = sa + c * (b - a)
    w = x
    v = w
    e = 0
    fx = f(x)
    fw = fx
    fv = fw

    iter = 0
    converged = false
    while (iter < max_iter && !converged)

        m = (sa + sb)/2
        t2 = 2 * tol

        if (abs(x - m) <= t2 - (sb - sa)/2)
            converged = true
        else
            r = 0
            q = r
            p = q

            if (tol < abs(e))

                r = (x - w) * (fx - fv)
                q = (x - v) * (fx - fw)
                p = (x - v) * q - (x - w) * r
                q = 2 * (q - r)

                if (0 < q)
                    p = - p
                end

                q = abs(q)

                r = e
                e = d

            end

            if (abs(p) < abs(q * r/2) && q * (sa - x) < p && p < q * (sb - x))
                d = p / q
                u = x + d

                if ((u - sa) < t2 || (sb - u) < t2)

                    if (x < m)
                        d = tol
                    else
                        d = - tol
                    end

                end

            else

                if (x < m)
                    e = sb - x
                else
                    e = sa - x
                end

                d = c * e

            end

            if (tol <= abs(d))
                u = x + d
            elseif (0 < d)
                u = x + tol
            else
                u = x - tol
            end

            fu = f(u)

            if (fu <= fx)

                if (u < x)
                    sb = x
                else
                    sa = x
                end

                v = w
                fv = fw
                w = x
                fw = fx
                x = u
                fx = fu

            else

                if (u < x)
                    sa = u
                else
                    sb = u
                end

                if (fu <= fw || w == x)
                    v = w
                    fv = fw
                    w = u
                    fw = fu
                elseif (fu <= fv || v == x || v == w)
                    v = u
                    fv = fu
                end

            end

        end

    iter += 1

    end

    return x
end
