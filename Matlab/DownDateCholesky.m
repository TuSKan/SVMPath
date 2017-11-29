   
function m_R = DownDateCholesky(m_R, T, idx)

   if (idx > 1)

        m_R(:,idx-1) = [];

        for i = idx - 1: size(m_R, 1) - 1,
            m_R(i:i+1,:) = givens(m_R(i, i), m_R(i+1, i)) * m_R(i:i+1,:);
        end

        m_R = m_R(1:end-1,:);

    else

        A = [-T(2) * T(3:end), 1; 
             eye(length(T(3:end))), zeros(length(T(3:end)),1)];

        m_R = m_R * A;

        for i = 1:size(m_R,1) - 1,
            m_R(i:i+1,:) = givens(m_R(i, i), m_R(i+1, i)) * m_R(i:i+1,:);
        end

        m_R = m_R(1:end-1,1:end-1);        

    end