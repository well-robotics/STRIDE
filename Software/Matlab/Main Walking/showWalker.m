function showWalker(q)
%SHOWRESULTS Summary of this function goes here
if size(q, 2)  == 1
    t = 0;
    anim = Animator.AmberAnimator(t, q(1:7));
    anim.pov = Animator.AnimatorPointOfView.West;
    anim.Animate(true);
    anim.isLooping = false;
    anim.updateWorldPosition = true;
    
else
    tvec = 1:size(q,2);
    q = q(1:7,:);
    anim = Animator.AmberAnimator(tvec, q);
    anim.pov = Animator.AnimatorPointOfView.West;
    anim.Animate(true);
    anim.isLooping = false;
    anim.updateWorldPosition = true;
    
    conGUI = Animator.AnimatorControls();
    conGUI.anim = anim;
end

