export interface ContextProviderInstance {
    constructor: typeof ContextProvider;
    config: ContextProviderConfig;
}
export interface ContextProviderConfig {
    name: ContextProviderName;
    domain?: string;
    entityName?: string;
    data?: object;
    extraData?: object;
    isEdited?: boolean;
    context?: object;
}
export declare enum ContextProviderName {
    PlanewaveCutoffDataManager = "PlanewaveCutoffDataManager",
    KGridFormDataManager = "KGridFormDataManager",
    QGridFormDataManager = "QGridFormDataManager",
    IGridFormDataManager = "IGridFormDataManager",
    QPathFormDataManager = "QPathFormDataManager",
    IPathFormDataManager = "IPathFormDataManager",
    KPathFormDataManager = "KPathFormDataManager",
    ExplicitKPathFormDataManager = "ExplicitKPathFormDataManager",
    ExplicitKPath2PIBAFormDataManager = "ExplicitKPath2PIBAFormDataManager",
    HubbardJContextManager = "HubbardJContextManager",
    HubbardUContextManager = "HubbardUContextManager",
    HubbardVContextManager = "HubbardVContextManager",
    HubbardContextManagerLegacy = "HubbardContextManagerLegacy",
    NEBFormDataManager = "NEBFormDataManager",
    BoundaryConditionsFormDataManager = "BoundaryConditionsFormDataManager",
    MLSettingsDataManager = "MLSettingsDataManager",
    MLTrainTestSplitDataManager = "MLTrainTestSplitDataManager",
    IonDynamicsContextProvider = "IonDynamicsContextProvider",
    CollinearMagnetizationDataManager = "CollinearMagnetizationDataManager",
    NonCollinearMagnetizationDataManager = "NonCollinearMagnetizationDataManager",
    QEPWXInputDataManager = "QEPWXInputDataManager",
    QENEBInputDataManager = "QENEBInputDataManager",
    VASPInputDataManager = "VASPInputDataManager",
    VASPNEBInputDataManager = "VASPNEBInputDataManager",
    NWChemInputDataManager = "NWChemInputDataManager"
}
export interface ContextProviderStatic {
    getConstructorConfig: (config: ContextProviderConfig) => ContextProviderInstance;
    createConfigFromContext: (config: ContextProviderConfig) => ContextProviderConfig;
    getExtraDataKeyByName: (name: string) => string;
    getIsEditedKeyByName: (name: string) => string;
}
export default class ContextProvider {
    config: ContextProviderConfig;
    name: ContextProviderName;
    domain?: string;
    entityName?: string;
    data?: object;
    extraData?: object;
    isEdited?: boolean;
    constructor(config: ContextProviderConfig);
    static getConstructorConfig(config: ContextProviderConfig): ContextProviderInstance;
    static createConfigFromContext(config: ContextProviderConfig): ContextProviderConfig & ({
        data: never;
        extraData: any;
        isEdited: any;
    } | {
        data?: undefined;
        extraData?: undefined;
        isEdited?: undefined;
    });
    setIsEdited(isEdited: boolean): void;
    getData(): void | object;
    setData(data: object): void;
    get defaultData(): void;
    transformData(data: object): object;
    yieldData(...transformDataArgs: any): {
        [x: string]: boolean | object | undefined;
    };
    yieldDataForRendering(): {
        [x: string]: boolean | object | undefined;
    };
    get extraDataKey(): string;
    static getExtraDataKeyByName(name: string): string;
    get isEditedKey(): string;
    static getIsEditedKeyByName(name: string): string;
    get isUnitContextProvider(): boolean;
    get isSubworkflowContextProvider(): boolean;
}
