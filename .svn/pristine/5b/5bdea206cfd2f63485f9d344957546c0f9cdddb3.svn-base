MODULE RivConnectivityMod
USE RivVarMod
IMPLICIT NONE

CONTAINS
    SUBROUTINE Label_clusters(selected_state, state, nclusters, cluster_labels)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: selected_state
        INTEGER, INTENT(IN) :: state(:,:)
        INTEGER, INTENT(OUT) :: nclusters
        INTEGER, INTENT(OUT) :: cluster_labels(:,:)
        INTEGER :: i,j
        
        DO i = 1,ns
            DO j = 1,nn
                cluster_labels(i,j) = 0
            ENDDO
        ENDDO
        
        nclusters = 0
        
        DO i = 1,ns
            DO j = 1,nn
                IF(Node_needs_label(selected_state, state(i,j), cluster_labels(i,j))) THEN
                    nclusters = nclusters+1
                    CALL Label_neighborhood(i,j,state(i,j), state, nclusters, cluster_labels)
                ENDIF
            ENDDO
        ENDDO

    END SUBROUTINE
    
    RECURSIVE SUBROUTINE Label_neighborhood(i, j, selected_state, state, label, cluster_labels)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: i, j, selected_state, state(:,:), label
        INTEGER, INTENT(OUT) :: cluster_labels(:,:)
!        Write(6,*) "i ", i, " j ", j
        cluster_labels(i,j) = label
!        IF(i.ne.1.and.i.ne.ns.and.j.ne.1.and.j.ne.nn) THEN
            IF(i.ne.1) THEN
        IF(Node_needs_label(selected_state, state(i-1, j), cluster_labels(i-1,j))) THEN
            CALL Label_neighborhood(i-1,j, selected_state, state, label, cluster_labels)
        ENDIF
            ENDIF
            if(i.ne.ns) then
        IF(Node_needs_label(selected_state, state(i+1, j), cluster_labels(i+1,j))) THEN
            CALL Label_neighborhood(i+1,j, selected_state, state, label, cluster_labels)
        ENDIF
            endif
            if(j.ne.1)then
        IF(Node_needs_label(selected_state, state(i, j-1), cluster_labels(i,j-1))) THEN
            CALL Label_neighborhood(i,j-1, selected_state, state, label, cluster_labels)
        ENDIF
             endif
   
            if(j.ne.nn)then
        IF(Node_needs_label(selected_state, state(i, j+1), cluster_labels(i,j+1))) THEN
            CALL Label_neighborhood(i,j+1, selected_state, state, label, cluster_labels)
        ENDIF
            endif
!        ENDIF
   
    END SUBROUTINE Label_neighborhood
    
    LOGICAL FUNCTION Node_needs_label(selected_state, node_state, node_label)
    IMPLICIT NONE
    INTEGER, INTENT (IN) :: selected_state, node_state, node_label
    
    IF(node_label == 0.and.((node_state.eq.selected_state.or.node_state.eq.4.or.node_state.eq.6).or.(selected_state ==-99))) THEN
        Node_needs_label = .TRUE.
    ELSE
        Node_needs_label = .FALSE.
    ENDIF
    
    END FUNCTION Node_needs_label
    
    SUBROUTINE Delete_clusters(ibc, cluster_labels, nclusters, nwetnodes)
    IMPLICIT NONE
    INTEGER, INTENT(INOUT) :: ibc(:,:), cluster_labels(:,:)
    INTEGER, INTENT(IN) :: nclusters
    INTEGER, INTENT(OUT) :: nwetnodes
    INTEGER :: i,j, status, max, maxindex
    INTEGER, ALLOCATABLE, DIMENSION(:) :: tmp
    ALLOCATE(tmp(nclusters), STAT = status)
    max = -1e6
    maxindex = 1
    tmp = 0
    nwetnodes = 0
    
        IF(nclusters > 1) THEN
            DO i = 1,ns
                DO j = 1,nn
                    IF(cluster_labels(i,j) > 0) THEN
                    tmp(cluster_labels(i,j)) = tmp(cluster_labels(i,j))+1
                    ENDIF
                ENDDO
            ENDDO
            DO i = 1,nclusters
                if(tmp(i) > max) THEN
                    max = tmp(i)
                    maxindex = i
                ENDIF
            ENDDO
            
            Write(6,*)'Delete Clusters'
            Write(6,*) nclusters
            Do i = 1,nclusters
                Write(6,*) 'Cluster ', i, ' Val ', tmp(i)
            ENDDO
            Write(6,*) maxindex
        ENDIF
        DO i = 1, ns
            DO j = 1,nn
                IF(nclusters > 1) THEN
                    IF(cluster_labels(i,j).gt.0.and.cluster_labels(i,j).ne.maxindex) THEN
                        ibc(i,j) = 0
                        hl(i,j) = 0.
			            u(i,j) = 0.
			            v(i,j) = 0.
                    ENDIF
                ENDIF
                If(ibc(i,j).ne.0) THEN
                    nwetnodes = nwetnodes+1
                ENDIF
            ENDDO
        ENDDO
        DEALLOCATE(tmp, STAT = status)
    END SUBROUTINE Delete_clusters

END MODULE