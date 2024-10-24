

createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        VonMises_YF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        VonMises_PF<
            BackStress<NullHardeningTensorFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        VonMises_YF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        VonMises_PF<
            BackStress<TensorLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        VonMises_YF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        VonMises_PF<
            BackStress<ArmstrongFrederickHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        VonMises_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        VonMises_PF<
            BackStress<NullHardeningTensorFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        VonMises_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        VonMises_PF<
            BackStress<TensorLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        VonMises_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        VonMises_PF<
            BackStress<ArmstrongFrederickHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        VonMises_YF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        DruckerPrager_PF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        VonMises_YF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        DruckerPrager_PF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        VonMises_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        DruckerPrager_PF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        VonMises_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        DruckerPrager_PF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        VonMises_YF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        ConstantDilatancy_PF<
            BackStress<NullHardeningTensorFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        VonMises_YF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        ConstantDilatancy_PF<
            BackStress<TensorLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        VonMises_YF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        ConstantDilatancy_PF<
            BackStress<ArmstrongFrederickHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        VonMises_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        ConstantDilatancy_PF<
            BackStress<NullHardeningTensorFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        VonMises_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        ConstantDilatancy_PF<
            BackStress<TensorLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        VonMises_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        ConstantDilatancy_PF<
            BackStress<ArmstrongFrederickHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        VonMises_YF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        MohrCoulomb_PF<
            BackStress<NullHardeningTensorFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        VonMises_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        MohrCoulomb_PF<
            BackStress<NullHardeningTensorFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        DruckerPrager_YF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        VonMises_PF<
            BackStress<NullHardeningTensorFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        DruckerPrager_YF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        VonMises_PF<
            BackStress<TensorLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        DruckerPrager_YF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        VonMises_PF<
            BackStress<ArmstrongFrederickHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        DruckerPrager_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        VonMises_PF<
            BackStress<NullHardeningTensorFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        DruckerPrager_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        VonMises_PF<
            BackStress<TensorLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        DruckerPrager_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        VonMises_PF<
            BackStress<ArmstrongFrederickHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        DruckerPrager_YF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        DruckerPrager_PF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        DruckerPrager_YF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        DruckerPrager_PF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        DruckerPrager_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        DruckerPrager_PF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        DruckerPrager_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        DruckerPrager_PF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        DruckerPrager_YF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        ConstantDilatancy_PF<
            BackStress<NullHardeningTensorFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        DruckerPrager_YF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        ConstantDilatancy_PF<
            BackStress<TensorLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        DruckerPrager_YF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        ConstantDilatancy_PF<
            BackStress<ArmstrongFrederickHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        DruckerPrager_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        ConstantDilatancy_PF<
            BackStress<NullHardeningTensorFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        DruckerPrager_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        ConstantDilatancy_PF<
            BackStress<TensorLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        DruckerPrager_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        ConstantDilatancy_PF<
            BackStress<ArmstrongFrederickHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        DruckerPrager_YF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        MohrCoulomb_PF<
            BackStress<NullHardeningTensorFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        DruckerPrager_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        MohrCoulomb_PF<
            BackStress<NullHardeningTensorFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        MohrCoulomb_YF<
            BackStress<NullHardeningTensorFunction>
            >, 
        VonMises_PF<
            BackStress<NullHardeningTensorFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        MohrCoulomb_YF<
            BackStress<NullHardeningTensorFunction>
            >, 
        VonMises_PF<
            BackStress<TensorLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        MohrCoulomb_YF<
            BackStress<NullHardeningTensorFunction>
            >, 
        VonMises_PF<
            BackStress<ArmstrongFrederickHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        MohrCoulomb_YF<
            BackStress<NullHardeningTensorFunction>
            >, 
        DruckerPrager_PF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        MohrCoulomb_YF<
            BackStress<NullHardeningTensorFunction>
            >, 
        DruckerPrager_PF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        MohrCoulomb_YF<
            BackStress<NullHardeningTensorFunction>
            >, 
        ConstantDilatancy_PF<
            BackStress<NullHardeningTensorFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        MohrCoulomb_YF<
            BackStress<NullHardeningTensorFunction>
            >, 
        ConstantDilatancy_PF<
            BackStress<TensorLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        MohrCoulomb_YF<
            BackStress<NullHardeningTensorFunction>
            >, 
        ConstantDilatancy_PF<
            BackStress<ArmstrongFrederickHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        MohrCoulomb_YF<
            BackStress<NullHardeningTensorFunction>
            >, 
        MohrCoulomb_PF<
            BackStress<NullHardeningTensorFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        TensionCutoff_YF<
            BackStress<NullHardeningTensorFunction>
            >, 
        VonMises_PF<
            BackStress<NullHardeningTensorFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        TensionCutoff_YF<
            BackStress<NullHardeningTensorFunction>
            >, 
        VonMises_PF<
            BackStress<TensorLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        TensionCutoff_YF<
            BackStress<NullHardeningTensorFunction>
            >, 
        VonMises_PF<
            BackStress<ArmstrongFrederickHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        TensionCutoff_YF<
            BackStress<NullHardeningTensorFunction>
            >, 
        DruckerPrager_PF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        TensionCutoff_YF<
            BackStress<NullHardeningTensorFunction>
            >, 
        DruckerPrager_PF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        TensionCutoff_YF<
            BackStress<NullHardeningTensorFunction>
            >, 
        ConstantDilatancy_PF<
            BackStress<NullHardeningTensorFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        TensionCutoff_YF<
            BackStress<NullHardeningTensorFunction>
            >, 
        ConstantDilatancy_PF<
            BackStress<TensorLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        TensionCutoff_YF<
            BackStress<NullHardeningTensorFunction>
            >, 
        ConstantDilatancy_PF<
            BackStress<ArmstrongFrederickHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        TensionCutoff_YF<
            BackStress<NullHardeningTensorFunction>
            >, 
        MohrCoulomb_PF<
            BackStress<NullHardeningTensorFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        DuncanChang_EL, 
        VonMises_YF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        VonMises_PF<
            BackStress<NullHardeningTensorFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        DuncanChang_EL, 
        VonMises_YF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        VonMises_PF<
            BackStress<TensorLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        DuncanChang_EL, 
        VonMises_YF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        VonMises_PF<
            BackStress<ArmstrongFrederickHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        DuncanChang_EL, 
        VonMises_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        VonMises_PF<
            BackStress<NullHardeningTensorFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        DuncanChang_EL, 
        VonMises_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        VonMises_PF<
            BackStress<TensorLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        DuncanChang_EL, 
        VonMises_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        VonMises_PF<
            BackStress<ArmstrongFrederickHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        DuncanChang_EL, 
        VonMises_YF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        DruckerPrager_PF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        DuncanChang_EL, 
        VonMises_YF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        DruckerPrager_PF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        DuncanChang_EL, 
        VonMises_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        DruckerPrager_PF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        DuncanChang_EL, 
        VonMises_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        DruckerPrager_PF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        DuncanChang_EL, 
        VonMises_YF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        ConstantDilatancy_PF<
            BackStress<NullHardeningTensorFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        DuncanChang_EL, 
        VonMises_YF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        ConstantDilatancy_PF<
            BackStress<TensorLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        DuncanChang_EL, 
        VonMises_YF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        ConstantDilatancy_PF<
            BackStress<ArmstrongFrederickHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        DuncanChang_EL, 
        VonMises_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        ConstantDilatancy_PF<
            BackStress<NullHardeningTensorFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        DuncanChang_EL, 
        VonMises_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        ConstantDilatancy_PF<
            BackStress<TensorLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        DuncanChang_EL, 
        VonMises_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        ConstantDilatancy_PF<
            BackStress<ArmstrongFrederickHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        DuncanChang_EL, 
        VonMises_YF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        MohrCoulomb_PF<
            BackStress<NullHardeningTensorFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        DuncanChang_EL, 
        VonMises_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        MohrCoulomb_PF<
            BackStress<NullHardeningTensorFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        DuncanChang_EL, 
        DruckerPrager_YF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        VonMises_PF<
            BackStress<NullHardeningTensorFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        DuncanChang_EL, 
        DruckerPrager_YF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        VonMises_PF<
            BackStress<TensorLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        DuncanChang_EL, 
        DruckerPrager_YF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        VonMises_PF<
            BackStress<ArmstrongFrederickHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        DuncanChang_EL, 
        DruckerPrager_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        VonMises_PF<
            BackStress<NullHardeningTensorFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        DuncanChang_EL, 
        DruckerPrager_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        VonMises_PF<
            BackStress<TensorLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        DuncanChang_EL, 
        DruckerPrager_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        VonMises_PF<
            BackStress<ArmstrongFrederickHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        DuncanChang_EL, 
        DruckerPrager_YF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        DruckerPrager_PF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        DuncanChang_EL, 
        DruckerPrager_YF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        DruckerPrager_PF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        DuncanChang_EL, 
        DruckerPrager_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        DruckerPrager_PF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        DuncanChang_EL, 
        DruckerPrager_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        DruckerPrager_PF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        DuncanChang_EL, 
        DruckerPrager_YF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        ConstantDilatancy_PF<
            BackStress<NullHardeningTensorFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        DuncanChang_EL, 
        DruckerPrager_YF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        ConstantDilatancy_PF<
            BackStress<TensorLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        DuncanChang_EL, 
        DruckerPrager_YF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        ConstantDilatancy_PF<
            BackStress<ArmstrongFrederickHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        DuncanChang_EL, 
        DruckerPrager_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        ConstantDilatancy_PF<
            BackStress<NullHardeningTensorFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        DuncanChang_EL, 
        DruckerPrager_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        ConstantDilatancy_PF<
            BackStress<TensorLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        DuncanChang_EL, 
        DruckerPrager_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        ConstantDilatancy_PF<
            BackStress<ArmstrongFrederickHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        DuncanChang_EL, 
        DruckerPrager_YF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        MohrCoulomb_PF<
            BackStress<NullHardeningTensorFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        DuncanChang_EL, 
        DruckerPrager_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        MohrCoulomb_PF<
            BackStress<NullHardeningTensorFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        DuncanChang_EL, 
        MohrCoulomb_YF<
            BackStress<NullHardeningTensorFunction>
            >, 
        VonMises_PF<
            BackStress<NullHardeningTensorFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        DuncanChang_EL, 
        MohrCoulomb_YF<
            BackStress<NullHardeningTensorFunction>
            >, 
        VonMises_PF<
            BackStress<TensorLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        DuncanChang_EL, 
        MohrCoulomb_YF<
            BackStress<NullHardeningTensorFunction>
            >, 
        VonMises_PF<
            BackStress<ArmstrongFrederickHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        DuncanChang_EL, 
        MohrCoulomb_YF<
            BackStress<NullHardeningTensorFunction>
            >, 
        DruckerPrager_PF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        DuncanChang_EL, 
        MohrCoulomb_YF<
            BackStress<NullHardeningTensorFunction>
            >, 
        DruckerPrager_PF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        DuncanChang_EL, 
        MohrCoulomb_YF<
            BackStress<NullHardeningTensorFunction>
            >, 
        ConstantDilatancy_PF<
            BackStress<NullHardeningTensorFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        DuncanChang_EL, 
        MohrCoulomb_YF<
            BackStress<NullHardeningTensorFunction>
            >, 
        ConstantDilatancy_PF<
            BackStress<TensorLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        DuncanChang_EL, 
        MohrCoulomb_YF<
            BackStress<NullHardeningTensorFunction>
            >, 
        ConstantDilatancy_PF<
            BackStress<ArmstrongFrederickHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        DuncanChang_EL, 
        MohrCoulomb_YF<
            BackStress<NullHardeningTensorFunction>
            >, 
        MohrCoulomb_PF<
            BackStress<NullHardeningTensorFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        DuncanChang_EL, 
        TensionCutoff_YF<
            BackStress<NullHardeningTensorFunction>
            >, 
        VonMises_PF<
            BackStress<NullHardeningTensorFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        DuncanChang_EL, 
        TensionCutoff_YF<
            BackStress<NullHardeningTensorFunction>
            >, 
        VonMises_PF<
            BackStress<TensorLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        DuncanChang_EL, 
        TensionCutoff_YF<
            BackStress<NullHardeningTensorFunction>
            >, 
        VonMises_PF<
            BackStress<ArmstrongFrederickHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        DuncanChang_EL, 
        TensionCutoff_YF<
            BackStress<NullHardeningTensorFunction>
            >, 
        DruckerPrager_PF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        DuncanChang_EL, 
        TensionCutoff_YF<
            BackStress<NullHardeningTensorFunction>
            >, 
        DruckerPrager_PF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        DuncanChang_EL, 
        TensionCutoff_YF<
            BackStress<NullHardeningTensorFunction>
            >, 
        ConstantDilatancy_PF<
            BackStress<NullHardeningTensorFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        DuncanChang_EL, 
        TensionCutoff_YF<
            BackStress<NullHardeningTensorFunction>
            >, 
        ConstantDilatancy_PF<
            BackStress<TensorLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        DuncanChang_EL, 
        TensionCutoff_YF<
            BackStress<NullHardeningTensorFunction>
            >, 
        ConstantDilatancy_PF<
            BackStress<ArmstrongFrederickHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        DuncanChang_EL, 
        TensionCutoff_YF<
            BackStress<NullHardeningTensorFunction>
            >, 
        MohrCoulomb_PF<
            BackStress<NullHardeningTensorFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);

